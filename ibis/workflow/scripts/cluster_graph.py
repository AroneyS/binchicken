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
    clusters = [elusive_edges.filter(pl.col("length") == 2).select("samples", "target_ids", sample = pl.col("samples"))]
    for n in range(3, MAX_COASSEMBLY_SAMPLES + 1):
        # Join pairs to n-1 clusters and add n clusters from elusive_edges
        clusters.append(
            clusters[-1]
            .explode("sample")
            .join(clusters[0].explode("sample"), on="sample", how="inner", suffix="_2")
            .select(
                pl.concat_list("samples", "samples_2").arr.unique(),
                pl.concat_list("target_ids", "target_ids_2").arr.unique(),
                n_way_edge = False,
                )
            .filter(pl.col("samples").arr.lengths() == n)
            .vstack(
                elusive_edges
                .filter(pl.col("length") == n)
                .select("samples", "target_ids", n_way_edge = True)
                )
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
        .sort("total_targets", descending=True)
        .with_columns(
            unique_samples = 
                pl.col("samples")
                .map(accumulate_clusters, return_dtype=pl.Boolean),
            )
        .filter(pl.col("unique_samples"))
    )

    if clusters.height == 0:
        return pl.DataFrame(schema=output_columns)

    # Find best recover samples
    import pdb; pdb.set_trace()




    clusters = clusters.with_columns(
        pl.col("subgraphview").apply(
            lambda x: x.edges(data=True)
            ).alias("edgeview"),
    ).select(
        pl.col("community").arr.sort().arr.join(",").alias("samples"),
        "length",
        pl.col("edgeview").apply(
            lambda x: sum([data["weight"] for _,_,data in x])
            ).alias("total_weight"),
        pl.col("edgeview").apply(
            lambda x: len(set([i for l in [data["target_ids"].split(",") for _,_,data in x] for i in l]))
            ).alias("total_targets"),
        "total_size",
        pl.col("community").apply(
            lambda x: umbrellas.filter(
                    pl.col("community").apply(
                        lambda y: all([(a in y) for a in x])
                    )
                ).head(1).get_column("recover_samples")
        ).flatten().alias("recover_samples"),
    )

    if clusters.height == 0:
        return pl.DataFrame(schema=output_columns)

    clusters = clusters.unique().with_columns(
        pl.col("samples").str.split(",").alias("samples_list"),
        pl.col("recover_samples").str.split(",").alias("recover_samples_list"),
    ).with_columns(
        pl.col("samples_list").apply(
            lambda x: elusive_edges.with_columns(
                    pl.col("sample1").is_in(x).alias("sample1_bool"),
                    pl.col("sample2").is_in(x).alias("sample2_bool"),
                ).filter(pl.col("sample1_bool") != pl.col("sample2_bool")
                ).with_columns(
                    pl.when(pl.col("sample1_bool")).then(pl.col("sample2")).otherwise(pl.col("sample1")).alias("other_sample"),
                )
            ).alias("relevant_edges"),
        pl.col("samples_list").apply(
            lambda x: elusive_edges.with_columns(
                    pl.col("sample1").is_in(x).alias("sample1_bool"),
                    pl.col("sample2").is_in(x).alias("sample2_bool"),
                ).filter(pl.col("sample1_bool") & pl.col("sample2_bool")
                ).select(
                    pl.col("target_ids").str.split(",").flatten().alias("target_ids"),
                ).get_column("target_ids")
            ).alias("relevant_targets"),
    ).with_columns(
        pl.when(
            (pl.col("recover_samples_list").arr.lengths() >= MAX_RECOVERY_SAMPLES) | (pl.col("relevant_edges").apply(lambda x: x.height) == 0)
        ).then(
            pl.Series("other_samples", [[]], dtype = pl.List(pl.Utf8))
        ).otherwise(
            pl.col("relevant_edges").apply(
                lambda x: x.with_columns(pl.col("target_ids").str.split(","))
                    .explode("target_ids")
                    .groupby("other_sample")
                    .agg(pl.count())
                    .sort(["count", "other_sample"], descending=True)
                    .get_column("other_sample")
                )
        ).alias("recover_candidates"),
    ).with_columns(
        pl.col("recover_samples_list")
            .arr.concat(pl.col("recover_candidates"))
            .alias("recover_samples"),
    )

    clusters = clusters.select([
        "samples", "length", "total_weight", "total_targets", "total_size", 
        pl.col("recover_samples")
            .apply(lambda s: s.to_frame("s").groupby("s", maintain_order=True).first().to_series())
            .arr.head(MAX_RECOVERY_SAMPLES)
            .arr.sort()
            .arr.join(",")
            .alias("recover_samples"),
    ]).sort(["total_targets", "samples"], descending=True).with_columns(
        pl.lit("coassembly_").alias("coassembly") + pl.arange(0, clusters.height).cast(pl.Utf8)
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
