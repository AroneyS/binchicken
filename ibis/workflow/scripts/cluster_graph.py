########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import polars as pl
import os
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms import components
from operator import itemgetter
import itertools

def pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    output_columns = ["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples", "coassembly"]

    if len(elusive_edges) == 0:
        return pl.DataFrame(schema=output_columns)

    elusive_edges = (
        elusive_edges
        .with_columns(
            pl.col("samples").str.split(",").arr.eval(pl.element().str.replace(r"\.1$", "")),
            pl.col("target_ids").str.split(","),
            )
        .with_columns(length = pl.col("samples").arr.lengths())
    )

    # Start with pair clusters
    clusters = [elusive_edges.lazy().filter(pl.col("length") == 2).select("samples", "target_ids", sample = pl.col("samples"))]
    for n in range(3, MAX_COASSEMBLY_SAMPLES + 1):
        # Join pairs to n-1 clusters and add n clusters from elusive_edges
        clusters.append(
            pl.concat([
                # n-1 way edges + matching edges
                clusters[-1]
                .explode("sample")
                .join(clusters[0].explode("sample"), on="sample", how="inner", suffix="_2")
                .select(
                    pl.concat_list("samples", "samples_2").arr.unique(),
                    pl.concat_list("target_ids", "target_ids_2").arr.unique(),
                    n_way_edge = False,
                    )
                .filter(pl.col("samples").arr.lengths() == n),
                # n-way edges
                elusive_edges
                .lazy()
                .filter(pl.col("length") == n)
                .select("samples", "target_ids", n_way_edge = True)
            ])
            .groupby(pl.col("samples").arr.sort().arr.join(","))
            .agg(
                pl.col("target_ids").flatten(),
                pl.any("n_way_edge"),
                )
            .filter(pl.col("n_way_edge"))
            .select(
                pl.col("samples").str.split(","),
                pl.col("target_ids").arr.unique(),
                sample = pl.col("samples").str.split(",")
                )
        )

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

        return pl.Series(choices)

    sample_targets = (
        elusive_edges
        .select("samples", "target_ids")
        .explode("samples")
        .explode("target_ids")
        .unique(["samples", "target_ids"])
    )

    # Find best clusters (each sample appearing only once per length)
    clusters = (
        pl.concat(clusters)
        .with_columns(
            pl.col("samples").arr.sort().arr.join(","),
            pl.col("target_ids").arr.sort().arr.join(","),
            total_targets = pl.col("target_ids").arr.lengths(),
            length = pl.col("samples").arr.lengths()
            )
        .explode("sample")
        .join(read_size.lazy(), on="sample", how="inner")
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
        .sort("total_targets", descending=True)
        .with_columns(
            unique_samples = 
                pl.col("samples")
                .map(accumulate_clusters, return_dtype=pl.Boolean),
            )
        .filter(pl.col("unique_samples"))
        .with_columns(
            recover_samples = pl.col("target_ids").str.split(",").apply(
                lambda x: sample_targets
                          .filter(pl.col("target_ids").is_in(x))
                          .groupby("samples")
                          .agg(pl.col("target_ids").unique().count())
                          .sort("target_ids", descending=True)
                          .head(MAX_RECOVERY_SAMPLES)
                          .get_column("samples")
            ).arr.sort().arr.join(","),
            )
        .with_row_count("coassembly")
        .select(
            "samples", "length", "total_targets", "total_size", "recover_samples",
            coassembly = pl.lit("coassembly_") + pl.col("coassembly").cast(pl.Utf8)
            )
        .collect()
    )

    return clusters

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

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
