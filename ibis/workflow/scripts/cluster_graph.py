########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import polars as pl
import itertools
import os
import logging

OUTPUT_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": int,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

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
        .select(
            pl.col("samples"),
            pl.col("samples_hash"),
            pl.col("length"),
            )
        .join(
            df2,
            on=["length"],
            how="left",
        )
        .filter(pl.col("samples").list.set_intersection("right_samples").list.lengths() >= pl.col("length"))
        .groupby("samples_hash")
        .agg(pl.col("extra_targets").flatten())
        .select(
            pl.col("samples_hash"),
            pl.col("extra_targets").list.unique().fill_null([]),
            )
    )

    return (
        df1
        .join(output, on="samples_hash", how="left")
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
        .join(samples_df, on="target_ids", how="left")
        .groupby("samples_hash", "recover_candidates")
        .agg(pl.count())
        .sort("count", descending=True)
        .groupby("samples_hash")
        .head(MAX_RECOVERY_SAMPLES)
        .groupby("samples_hash")
        .agg("recover_candidates")
    )

    return df.join(output, on="samples_hash", how="left")

def pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20,
        MIN_CLUSTER_TARGETS=1,
        MAX_SAMPLES_COMBINATIONS=100,
        EXCLUDE_COASSEMBLIES=[]):

    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    if len(elusive_edges) == 0:
        logging.warning("No elusive edges found")
        return pl.DataFrame(schema=OUTPUT_COLUMNS)

    is_pooled = any(elusive_edges["style"] == "pool")

    with pl.StringCache():
        elusive_edges = (
            elusive_edges
            .lazy()
            .with_columns(
                pl.col("samples")
                    .str.split(",")
                    .list.eval(pl.element().str.replace(r"\.1$", ""))
                    .cast(pl.List(pl.Categorical)),
                pl.col("target_ids")
                    .str.split(",")
                    .cast(pl.List(pl.UInt32)),
                )
            .with_columns(
                samples_hash = pl.col("samples").list.sort().hash(),
                length = pl.col("samples").list.lengths()
                )
        )

        read_size = (
            read_size
            .lazy()
            .select(
                "read_size",
                samples = pl.col("sample").cast(pl.Categorical),
                )
        )

        if EXCLUDE_COASSEMBLIES:
            excluded_coassemblies = (
                pl.LazyFrame({"samples": EXCLUDE_COASSEMBLIES})
                .with_columns(
                    samples_hash = pl.col("samples")
                        .str.split(",")
                        .cast(pl.List(pl.Categorical))
                        .list.sort().hash()
                    )
            )
        else:
            excluded_coassemblies = pl.LazyFrame(schema={"samples_hash": pl.UInt64})

        if MAX_COASSEMBLY_SAMPLES == 1:
            logging.info("Skipping clustering, using single-sample clusters")
            clusters = [
                elusive_edges
                .explode("samples")
                .groupby("samples")
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
                elusive_edges
                .filter(pl.col("style") == "match")
                .filter(pl.col("cluster_size") >= MIN_COASSEMBLY_SAMPLES)
                .filter(pl.col("target_ids").list.lengths() >= MIN_CLUSTER_TARGETS)
                .select("samples", "target_ids", samples_hash = pl.col("samples").list.sort().hash())
            ]

            if is_pooled:
                clusters.append(
                    elusive_edges
                    .filter(pl.col("style") == "pool")
                    .filter(pl.col("cluster_size") >= MIN_COASSEMBLY_SAMPLES)
                    # Prevent combinatorial explosion (also, large clusters are less useful for distinguishing between clusters)
                    .filter(pl.col("samples").list.lengths() < MAX_SAMPLES_COMBINATIONS)
                    .with_columns(
                        samples_combinations = pl.struct(["length", "cluster_size"]).apply(
                            lambda x: [i for i in itertools.combinations(range(x["length"]), x["cluster_size"])],
                            return_dtype=pl.List(pl.List(pl.Int64)),
                            )
                        )
                    .explode("samples_combinations")
                    .with_columns(pl.col("samples").list.take(pl.col("samples_combinations")))
                    .with_columns(samples_hash = pl.col("samples").list.sort().hash())
                    .select("samples_hash", "samples", "target_ids", "cluster_size")
                    .groupby("samples_hash")
                    .agg(
                        pl.first("samples"),
                        pl.col("target_ids").flatten(),
                        pl.first("cluster_size"),
                        )
                    .filter(pl.col("target_ids").list.lengths() >= MIN_CLUSTER_TARGETS)
                    .select("samples", "target_ids", "samples_hash")
                )

        sample_targets = (
            elusive_edges
            .select("target_ids", recover_candidates = pl.col("samples"))
            .explode("recover_candidates")
            .explode("target_ids")
            .unique()
            .groupby("recover_candidates")
            .agg("target_ids")
        )

        logging.info("Filtering clusters (each sample restricted to only once per cluster size)")
        clusters = (
            pl.concat(clusters)
            .join(excluded_coassemblies, on="samples_hash", how="anti")
            .with_columns(length = pl.col("samples").list.lengths())
            .explode("samples")
            .join(read_size, on="samples", how="left")
            .groupby("samples_hash")
            .agg(
                pl.col("samples"),
                pl.first("length"),
                pl.first("target_ids"),
                total_size = pl.sum("read_size"),
                )
            .filter(
                pl.when(MAX_COASSEMBLY_SIZE is not None)
                .then(pl.col("total_size") <= MAX_COASSEMBLY_SIZE)
                .otherwise(True)
                )
            .with_columns(
                total_targets = pl.col("target_ids").list.lengths(),
            )
            .sort("total_targets", "total_size", descending=[True, False])
            .with_columns(
                unique_samples = 
                    pl.col("samples")
                    .map(accumulate_clusters, return_dtype=pl.Boolean),
                )
            .filter(pl.col("unique_samples"))
            .drop("unique_samples")
            .collect(streaming=True)
            .pipe(
                join_list_subsets,
                df2=elusive_edges
                    .filter(pl.col("style") == "pool")
                    .filter(pl.col("samples").list.lengths() >= MAX_SAMPLES_COMBINATIONS)
                    .collect(streaming=True),
                )
            .with_columns(pl.concat_list("target_ids", "extra_targets").list.unique())
            .pipe(
                find_recover_candidates,
                sample_targets.collect(streaming=True),
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
                total_targets = pl.col("target_ids").list.lengths(),
                )
            .sort("total_targets", "total_size", descending=[True, False])
            .with_row_count("coassembly")
            .select(
                "samples", "length", "total_targets", "total_size", "recover_samples",
                coassembly = pl.lit("coassembly_") + pl.col("coassembly").cast(pl.Utf8)
                )
        )

    logging.info("Done")

    return clusters

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    MAX_COASSEMBLY_SIZE = snakemake.params.max_coassembly_size * 10**9 if snakemake.params.max_coassembly_size else None
    MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
    MIN_COASSEMBLY_SAMPLES = snakemake.params.num_coassembly_samples
    MAX_RECOVERY_SAMPLES = snakemake.params.max_recovery_samples
    EXCLUDE_COASSEMBLIES = snakemake.params.exclude_coassemblies
    elusive_edges_path = snakemake.input.elusive_edges
    read_size_path = snakemake.input.read_size
    elusive_clusters_path = snakemake.output.elusive_clusters

    elusive_edges = pl.read_csv(elusive_edges_path, separator="\t", dtypes={"target_ids": str})
    read_size = pl.read_csv(read_size_path, has_header=False, new_columns=["sample", "read_size"])

    if elusive_edges.height > 10**4:
        min_cluster_targets = 10
    else:
        min_cluster_targets = 1

    clusters = pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES,
        EXCLUDE_COASSEMBLIES=EXCLUDE_COASSEMBLIES,
        MIN_CLUSTER_TARGETS=min_cluster_targets,
        MAX_SAMPLES_COMBINATIONS=100,
        )
    clusters.write_csv(elusive_clusters_path, separator="\t")
