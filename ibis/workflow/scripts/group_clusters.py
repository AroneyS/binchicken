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
            edge.with_columns(precluster = pl.lit(precluster))
        )

    return (
        pl.concat(dfs)
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

def targets_processing(targets, preclusters):
    if len(targets) == 1:
        return targets[0]

    dfs = []
    for target, precluster in zip(targets, preclusters):
        dfs.append(
            target.with_columns(precluster = pl.lit(precluster))
        )

    return (
        pl.concat(dfs)
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

def cluster_processing(clusters, preclusters):
    if len(clusters) == 1:
        return clusters[0]

    return pl.concat(clusters)

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl
    print(f"Polars using {str(pl.threadpool_size())} threads")

    input_edges = snakemake.input.elusive_edges
    input_targets = snakemake.input.targets
    input_clusters = snakemake.input.elusive_clusters
    output_edges_path = snakemake.output.elusive_edges
    output_targets_path = snakemake.output.targets
    output_clusters_path = snakemake.output.elusive_clusters
    preclusters = snakemake.params.preclusters

    edges = []
    for edge in input_edges:
        edges.append(pl.scan_csv(edge, separator="\t", dtypes=EDGES_COLUMNS))
    output_edges = edges_processing(edges, preclusters)
    output_edges.collect().write_csv(output_edges_path, separator="\t")

    targets = []
    for target in input_targets:
        targets.append(pl.scan_csv(target, separator="\t", dtypes=TARGETS_COLUMNS))
    output_targets = targets_processing(targets, preclusters)
    output_targets.collect().write_csv(output_targets_path, separator="\t")

    clusters = []
    for cluster in input_clusters:
        clusters.append(pl.scan_csv(cluster, separator="\t", dtypes=ELUSIVE_CLUSTERS_COLUMNS))
    output_clusters = cluster_processing(clusters, preclusters)
    output_clusters.collect().write_csv(output_clusters_path, separator="\t")
