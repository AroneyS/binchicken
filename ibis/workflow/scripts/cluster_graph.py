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

def pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20):

    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    if len(elusive_edges) == 0:
        logging.warning("No elusive edges found")
        return pl.DataFrame(schema=OUTPUT_COLUMNS)

    elusive_edges = (
        elusive_edges
        .with_columns(
            pl.col("samples").str.split(",").list.eval(pl.element().str.replace(r"\.1$", "")),
            pl.col("target_ids").str.split(","),
            )
        .with_columns(length = pl.col("samples").list.lengths())
    )

    if MAX_COASSEMBLY_SAMPLES == 1:
        logging.info("Skipping clustering, using single-sample clusters")
        clusters = (
            elusive_edges
            .explode("samples")
            .groupby("samples")
            .agg(pl.col("target_ids").flatten())
            .with_columns(
                pl.col("samples").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)),
                pl.col("target_ids").list.sort().list.unique(),
                sample = pl.col("samples").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)),
            )
        )
    else:
        logging.info("Forming candidate sample clusters")
        clusters = pl.concat([
            # Pair clusters
            elusive_edges
            .filter(pl.col("style") == "match")
            .filter(pl.col("length") == 2)
            .select("samples", "target_ids", sample = pl.col("samples")),
            # 3+ way clusters
            elusive_edges
            .filter(pl.col("style") == "pool")
            .with_columns(
                samples = pl.struct(["samples", "cluster_size"]).apply(
                    lambda x: [i for i in itertools.combinations(x["samples"], x["cluster_size"])],
                    return_dtype=pl.List(pl.List(pl.Utf8)),
                    )
                )
            .explode("samples")
            .select(pl.col("samples").list.sort().list.join(","), "target_ids")
            .groupby("samples")
            .agg(pl.col("target_ids").flatten())
            .with_columns(
                pl.col("samples").str.split(",")
            )
            .with_columns(
                two_way_targets = pl.col("samples").apply(
                    lambda x: elusive_edges
                            .filter(pl.col("style") == "match")
                            .filter(pl.col("samples").list.eval(pl.element().is_in(x)).list.sum() == 2)
                            .select(pl.col("target_ids").flatten())
                            .get_column("target_ids")
                            .to_list(),
                    skip_nulls=False
                    ).cast(pl.List(pl.Utf8)),
                )
            .select(
                "samples",
                pl.concat_list("target_ids", "two_way_targets").list.sort().list.unique(),
                sample = pl.col("samples"),
                )
        ])

    def accumulate_clusters(x):
        clustered_samples = []
        choices = []
        for cluster in x:
            cluster_samples = cluster.split(",")
            if all([x not in clustered_samples for x in cluster_samples]):
                choices.append(True)
                clustered_samples += cluster_samples
            else:
                choices.append(False)

        return pl.Series(choices, dtype=pl.Boolean)

    sample_targets = (
        elusive_edges
        .select("samples", "target_ids")
        .explode("samples")
        .explode("target_ids")
        .unique(["samples", "target_ids"])
    )

    logging.info("Filtering clusters (each sample restricted to only once per cluster size)")
    clusters = (
        clusters
        .with_columns(
            pl.col("samples").list.sort().list.join(","),
            pl.col("target_ids").list.sort().list.join(","),
            total_targets = pl.col("target_ids").list.lengths(),
            length = pl.col("samples").list.lengths()
            )
        .explode("sample")
        .join(read_size, on="sample", how="inner")
        .groupby("samples")
        .agg(
            pl.first("length"),
            pl.first("target_ids"),
            pl.first("total_targets"),
            total_size = pl.sum("read_size"),
            )
        .filter(pl.col("length") >= MIN_COASSEMBLY_SAMPLES)
        .filter(
            pl.when(MAX_COASSEMBLY_SIZE is not None)
            .then(pl.col("total_size") <= MAX_COASSEMBLY_SIZE)
            .otherwise(True)
            )
        .sort("total_targets", "samples", descending=True)
        .with_columns(
            unique_samples = 
                pl.col("samples")
                .map(accumulate_clusters, return_dtype=pl.Boolean),
            )
        .filter(pl.col("unique_samples"))
        .with_columns(
            recover_candidates = pl.col("target_ids").str.split(",").apply(
                lambda x: sample_targets
                          .filter(pl.col("target_ids").is_in(x))
                          .groupby("samples")
                          .agg(pl.col("target_ids").unique().count())
                          .sort("target_ids", descending=True)
                          .head(MAX_RECOVERY_SAMPLES)
                          .get_column("samples"),
                return_dtype=pl.List(pl.Utf8),
            ),
            )
        .with_columns(
            recover_samples = pl.concat_list(pl.col("samples").str.split(","), "recover_candidates")
            .list.unique(maintain_order=True)
            .list.head(MAX_RECOVERY_SAMPLES)
            .list.sort()
            .list.join(","),
            )
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
    elusive_edges_path = snakemake.input.elusive_edges
    read_size_path = snakemake.input.read_size
    elusive_clusters_path = snakemake.output.elusive_clusters

    elusive_edges = pl.read_csv(elusive_edges_path, separator="\t", dtypes={"target_ids": str})
    read_size = pl.read_csv(read_size_path, has_header=False, new_columns=["sample", "read_size"])

    clusters = pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES
        )
    clusters.write_csv(elusive_clusters_path, separator="\t")
