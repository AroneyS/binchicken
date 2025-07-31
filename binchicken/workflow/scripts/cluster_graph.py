#!/usr/bin/env python3
########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import polars as pl
import itertools
import os
import logging
import argparse

OUTPUT_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

TARGET_WEIGHTING_COLUMNS = {
    "target": str,
    "weight": float,
}

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/cluster_graph.py \
#   --elusive-edges {input.elusive_edges} \
#   --read-size {input.read_size} \
#   --targets-weighted {input.targets_weighted} \
#   --elusive-clusters {output.elusive_clusters} \
#   --max-coassembly-size {params.max_coassembly_size} \
#   --max-coassembly-samples {params.max_coassembly_samples} \
#   --num-coassembly-samples {params.num_coassembly_samples} \
#   --max-recovery-samples {params.max_recovery_samples} \
#   --max-samples-combinations {params.max_samples_combinations} \
#   --coassembly-samples {input.coassembly_samples} \
#   --exclude-coassemblies {input.exclude_coassemblies} \
#   --single-assembly {params.single_assembly} \
#   --anchor-samples {input.anchor_samples} \
#   --threads {threads} \
#   --log {log}
# """

def join_list_subsets(df1, df2):
    """
    Join two DataFrames/LazyFrames by strict subset of list column
    """
    df2 = (
        df2.select(
            right_samples = pl.col("samples"),
            extra_targets = pl.col("target_ids"),
            length = pl.col("cluster_size").cast(pl.UInt32),
            )
    )

    output = (
        df1
        .select("samples", "samples_hash", "length")
        .join(df2, on="length")
        .filter(pl.col("right_samples").is_not_null())
        .filter(pl.col("samples").list.set_intersection("right_samples").list.len() >= pl.col("length"))
        .group_by("samples_hash")
        .agg(pl.col("extra_targets").flatten())
        .select(
            pl.col("samples_hash"),
            pl.col("extra_targets").list.unique().fill_null([]),
            )
    )

    return (
        df1
        .join(output, on="samples_hash", how="left", coalesce=True)
        .with_columns(pl.col("extra_targets").fill_null([]))
        )

def accumulate_clusters(x):
    clustered_samples = {}
    choices = []
    for cluster in x:
        cluster_size = len(cluster)
        if cluster_size not in clustered_samples:
            clustered_samples[cluster_size] = set()
        if all([s not in clustered_samples[cluster_size] for s in cluster]):
            choices.append(True)
            clustered_samples[cluster_size].update(cluster)
        else:
            choices.append(False)

    return pl.Series(choices, dtype=pl.Boolean)

def find_recover_candidates(df, samples_df, MAX_RECOVERY_SAMPLES=20):
    samples_df = samples_df.explode("target_ids")

    output = (
        df
        .select("samples_hash", "target_ids")
        .explode("target_ids")
        .join(samples_df, on="target_ids", how="left", coalesce=True)
        .group_by("samples_hash", "recover_candidates")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
        .group_by("samples_hash")
        .head(MAX_RECOVERY_SAMPLES)
        .group_by("samples_hash")
        .agg("recover_candidates")
    )

    return df.join(output, on="samples_hash", how="left", coalesce=True)

def pipeline(
        elusive_edges,
        read_size,
        weightings=None,
        anchor_samples=set(),
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20,
        MIN_CLUSTER_TARGETS=1,
        MAX_SAMPLES_COMBINATIONS=100,
        COASSEMBLY_SAMPLES=[],
        EXCLUDE_COASSEMBLIES=[],
        single_assembly=False):

    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(elusive_edges) == 0:
        logging.warning("No elusive edges found")
        return pl.DataFrame(orient="row", schema=OUTPUT_COLUMNS)

    is_pooled = any(elusive_edges["style"] == "pool")

    with pl.StringCache():
        elusive_edges = (
            elusive_edges
            .with_columns(
                pl.col("samples")
                    .str.split(",")
                    .cast(pl.List(pl.Categorical)),
                pl.col("target_ids")
                    .str.split(",")
                    .cast(pl.List(pl.UInt64)),
                )
            .with_columns(
                samples_hash = pl.col("samples").list.sort().hash(),
                length = pl.col("samples").list.len()
                )
        )

        if weightings is not None:
            if weightings.height == 0:
                logging.error("No target weightings found")
                return pl.DataFrame(schema=OUTPUT_COLUMNS)

            weightings = (
                weightings
                .select(target_ids = pl.col("target").cast(pl.UInt64), weight = "weight")
            )
            weightings_dict = dict(weightings.iter_rows())
        else:
            weightings_dict = {}

        if COASSEMBLY_SAMPLES:
            coassembly_edges = (
                elusive_edges
                .with_columns(
                    pl.col("samples").list.eval(pl.element().filter(pl.element().is_in(COASSEMBLY_SAMPLES)))
                    )
                .with_columns(length = pl.col("samples").list.len())
                .filter(pl.col("samples").list.len() >= pl.col("cluster_size"))
            )
        else:
            coassembly_edges = elusive_edges

        read_size = (
            read_size
            .select(
                "read_size",
                samples = pl.col("sample").cast(pl.Categorical),
                )
        )

        if EXCLUDE_COASSEMBLIES:
            excluded_coassemblies = (
                pl.DataFrame({"samples": EXCLUDE_COASSEMBLIES})
                .with_columns(
                    samples_hash = pl.col("samples")
                        .str.split(",")
                        .cast(pl.List(pl.Categorical))
                        .list.sort().hash()
                    )
            )
        else:
            excluded_coassemblies = pl.DataFrame(schema={"samples_hash": pl.UInt64})

        if MAX_COASSEMBLY_SAMPLES == 1:
            logging.info("Skipping clustering, using single-sample clusters")
            clusters = [
                elusive_edges
                .explode("samples")
                .filter((not COASSEMBLY_SAMPLES) | pl.col("samples").is_in(COASSEMBLY_SAMPLES))
                .group_by("samples")
                .agg(pl.col("target_ids").flatten())
                .with_columns(
                    pl.concat_list(pl.col("samples")),
                    pl.col("target_ids").list.sort().list.unique(),
                    samples_hash = pl.concat_list(pl.col("samples")).list.sort().hash(),
                )
            ]
        else:
            logging.info("Forming candidate sample clusters")
            clusters = [
                coassembly_edges
                .filter(pl.col("style") == "match")
                .filter(pl.col("cluster_size") >= MIN_COASSEMBLY_SAMPLES)
                .filter(pl.col("target_ids").list.len() >= MIN_CLUSTER_TARGETS)
                .select("samples", "target_ids", samples_hash = pl.col("samples").list.sort().hash())
            ]

            if is_pooled:
                clusters.append(
                    coassembly_edges
                    .filter(pl.col("style") == "pool")
                    .filter(pl.col("cluster_size") >= MIN_COASSEMBLY_SAMPLES)
                    # Prevent combinatorial explosion (also, large clusters are less useful for distinguishing between clusters)
                    .filter(pl.col("samples").list.len() <= MAX_SAMPLES_COMBINATIONS)
                    .with_columns(
                        samples_combinations = pl.struct(["length", "cluster_size"]).map_elements(
                            lambda x: [i for i in itertools.combinations(range(x["length"]), x["cluster_size"])],
                            return_dtype=pl.List(pl.List(pl.Int64)),
                            )
                        )
                    .explode("samples_combinations")
                    .with_columns(pl.col("samples").list.gather(pl.col("samples_combinations")))
                    .with_columns(samples_hash = pl.col("samples").list.sort().hash())
                    .select("samples_hash", "samples", "target_ids", "cluster_size")
                    .group_by("samples_hash")
                    .agg(
                        pl.first("samples"),
                        pl.col("target_ids").flatten(),
                        pl.first("cluster_size"),
                        )
                    .filter(pl.col("target_ids").list.len() >= MIN_CLUSTER_TARGETS)
                    .select("samples", "target_ids", "samples_hash")
                )

        sample_targets = (
            elusive_edges
            .select("target_ids", recover_candidates = pl.col("samples"))
            .explode("recover_candidates")
            .explode("target_ids")
            .unique()
            .group_by("recover_candidates")
            .agg("target_ids")
        )

        def filter_max_coassembly_size(df, MAX_COASSEMBLY_SIZE):
            if MAX_COASSEMBLY_SIZE is None:
                return df
            else:
                return df.filter(pl.col("total_size") <= MAX_COASSEMBLY_SIZE)

        logging.info("Filtering clusters (each sample restricted to only once per cluster size)")
        clusters = (
            pl.concat(clusters)
            .join(excluded_coassemblies, on="samples_hash", how="anti")
            .filter((not anchor_samples) | (pl.col("samples").list.eval(pl.element().is_in(anchor_samples)).list.any()))
            .with_columns(length = pl.col("samples").list.len())
            .explode("samples")
            .join(read_size, on="samples", how="left", coalesce=True)
            .group_by("samples_hash")
            .agg(
                pl.col("samples"),
                pl.first("length"),
                pl.first("target_ids"),
                total_size = pl.sum("read_size"),
                )
            .pipe(
                filter_max_coassembly_size,
                MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
                )
            .with_columns(weighting = pl.lit(weightings is not None))
            .with_columns(
                total_targets = pl.when(pl.col("weighting"))
                .then(pl.col("target_ids").list.eval(pl.element().replace_strict(weightings_dict, default=0)).list.sum())
                .otherwise(pl.col("target_ids").list.len()),
            )
            .sort("total_targets", "total_size", descending=[True, False])
            .with_columns(
                unique_samples = 
                    pl.col("samples")
                    .map_batches(accumulate_clusters, return_dtype=pl.Boolean),
                )
            .filter(pl.col("unique_samples"))
            .drop("unique_samples")
            .pipe(
                join_list_subsets,
                df2=coassembly_edges
                    .filter(pl.col("style") == "pool")
                    .filter(pl.col("samples").list.len() > MAX_SAMPLES_COMBINATIONS),
                )
            .with_columns(pl.concat_list("target_ids", "extra_targets").list.unique())
            .pipe(
                find_recover_candidates,
                sample_targets,
                MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES,
                )
            .with_columns(
                pl.col("samples").cast(pl.List(pl.Utf8)).list.sort().list.join(","),
                recover_samples = pl.concat_list(pl.col("samples"), "recover_candidates")
                .list.unique(maintain_order=True)
                .list.head(MAX_RECOVERY_SAMPLES)
                .cast(pl.List(pl.Utf8))
                .list.sort()
                .list.join(","),
                total_targets = pl.when(pl.col("weighting"))
                    .then(pl.col("total_targets"))
                    .otherwise(pl.col("target_ids").list.len()),
                )
            .sort("total_targets", "total_size", descending=[True, False])
            .with_row_index("coassembly")
            .select(
                "samples", "length", "total_targets", "total_size", "recover_samples",
                coassembly = pl.when(single_assembly)
                    .then(pl.col("samples"))
                    .otherwise(pl.lit("coassembly_") + pl.col("coassembly").cast(pl.Utf8))
                )
        )

    logging.info("Done")

    return clusters

def main():
    parser = argparse.ArgumentParser(description="Cluster graph pipeline script.")
    parser.add_argument("--elusive-edges", required=True, help="Path to elusive edges input file")
    parser.add_argument("--read-size", required=True, help="Path to read size input file")
    parser.add_argument("--targets-weighted", nargs='?', const=None, default=None, help="Path to targets weighted input file (optional)")
    parser.add_argument("--elusive-clusters", required=True, help="Path to output elusive clusters file")
    parser.add_argument("--max-coassembly-size", type=float, nargs='?', const=None, default=None, help="Max coassembly size in GB (optional)")
    parser.add_argument("--max-coassembly-samples", type=int, default=2, help="Max coassembly samples")
    parser.add_argument("--num-coassembly-samples", type=int, default=2, help="Num coassembly samples (min)")
    parser.add_argument("--max-recovery-samples", type=int, default=20, help="Max recovery samples")
    parser.add_argument("--max-samples-combinations", type=int, default=100, help="Max samples combinations for pooled clusters")
    parser.add_argument("--coassembly-samples", nargs='?', default=None, help="List file of coassembly samples")
    parser.add_argument("--exclude-coassemblies", nargs='?', default=None, help="List file of excluded coassemblies")
    parser.add_argument("--single-assembly", type=lambda x: x.lower() == 'true', nargs='?', const=True, default=False, help="Single assembly mode (True/False)")
    parser.add_argument("--anchor-samples", nargs='?', default=None, help="List file of anchor samples")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if args.log:
        logging.basicConfig(
            filename=args.log,
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    MAX_COASSEMBLY_SIZE = args.max_coassembly_size * 10**9 if args.max_coassembly_size else None
    MAX_COASSEMBLY_SAMPLES = args.max_coassembly_samples
    MIN_COASSEMBLY_SAMPLES = args.num_coassembly_samples
    MAX_RECOVERY_SAMPLES = args.max_recovery_samples
    MAX_SAMPLES_COMBINATIONS = args.max_samples_combinations
    COASSEMBLY_SAMPLES = pl.read_csv(args.coassembly_samples, has_header=False, new_columns=["sample"]).get_column("sample").to_list() if args.coassembly_samples else []
    EXCLUDE_COASSEMBLIES = pl.read_csv(args.exclude_coassemblies, separator="\t", has_header=False, new_columns=["coassembly"]).get_column("coassembly").to_list() if args.exclude_coassemblies else []
    single_assembly = args.single_assembly
    anchor_samples = set(pl.read_csv(args.anchor_samples, has_header=False, new_columns=["sample"]).get_column("sample").to_list()) if args.anchor_samples else set()

    elusive_edges = pl.read_csv(args.elusive_edges, separator="\t", schema_overrides={"target_ids": str})
    read_size = pl.read_csv(args.read_size, has_header=False, new_columns=["sample", "read_size"])

    if args.targets_weighted:
        weightings = pl.read_csv(args.targets_weighted, separator="\t", schema_overrides=TARGET_WEIGHTING_COLUMNS)
    else:
        weightings = None

    if elusive_edges.height > 10**4:
        min_cluster_targets = 10
    else:
        min_cluster_targets = 1

    clusters = pipeline(
        elusive_edges,
        read_size,
        weightings,
        anchor_samples=anchor_samples,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES,
        COASSEMBLY_SAMPLES=COASSEMBLY_SAMPLES,
        EXCLUDE_COASSEMBLIES=EXCLUDE_COASSEMBLIES,
        MIN_CLUSTER_TARGETS=min_cluster_targets,
        MAX_SAMPLES_COMBINATIONS=MAX_SAMPLES_COMBINATIONS,
        single_assembly=single_assembly,
    )
    clusters.write_csv(args.elusive_clusters, separator="\t")

if __name__ == "__main__":
    main()
