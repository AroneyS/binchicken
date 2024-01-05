#########################
### group_clusters.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os

EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
    }

TARGETS_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "target": str,
    }

ELUSIVE_CLUSTERS_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": int,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

def edges_processing(edges, preclusters):
    if len(edges) == 1:
        return edges[0]

    dfs = []
    for edge, precluster in zip(edges, preclusters):
        dfs.append(
            edge
            .with_columns(precluster = pl.lit(precluster))
            .with_columns(pl.col("target_ids").str.split(","))
            .explode("target_ids")
            .select(
                "style", "cluster_size", "samples",
                target_ids = pl.concat_str("precluster", pl.lit("_"), "target_ids")
                )
            .group_by("style", "cluster_size", "samples")
            .agg(pl.col("target_ids"))
            .with_columns(pl.col("target_ids").list.join(","))
        )

    return (
        pl.concat(dfs)
    )

def targets_processing(targets, preclusters):
    if len(targets) == 1:
        return targets[0]

    dfs = []
    for target, precluster in zip(targets, preclusters):
        dfs.append(
            target
            .with_columns(precluster = pl.lit(precluster))
            .with_columns(pl.col("target").str.split(","))
            .explode("target")
            .select(
                "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
                target = pl.concat_str("precluster", pl.lit("_"), "target")
                )
            .group_by("gene", "sample", "sequence", "num_hits", "coverage", "taxonomy")
            .agg(pl.col("target"))
            .with_columns(pl.col("target").list.join(","))
        )

    return (
        pl.concat(dfs)
    )

def cluster_processing(clusters, _):
    if len(clusters) == 1:
        return clusters[0]

    return (
        pl.concat(clusters)
        .drop("coassembly")
        .sort("total_targets", "total_size", descending=[True, False])
        .with_row_count("coassembly")
        .select(
            "samples", "length", "total_targets", "total_size", "recover_samples",
            coassembly = pl.lit("coassembly_") + pl.col("coassembly").cast(pl.Utf8)
            )
    )

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl
    print(f"Polars using {str(pl.threadpool_size())} threads")

    input_clusters = snakemake.input.elusive_clusters
    output_clusters_path = snakemake.output.elusive_clusters
    preclusters = snakemake.params.preclusters

    clusters = []
    for cluster in input_clusters:
        clusters.append(pl.scan_csv(cluster, separator="\t", dtypes=ELUSIVE_CLUSTERS_COLUMNS))
    output_clusters = cluster_processing(clusters, preclusters)
    output_clusters.collect().write_csv(output_clusters_path, separator="\t")
