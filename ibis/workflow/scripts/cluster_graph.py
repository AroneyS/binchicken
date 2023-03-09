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

    elusive_edges = elusive_edges.with_columns(
        pl.col("sample1").str.replace(r"\.1$", ""),
        pl.col("sample2").str.replace(r"\.1$", ""),
        )

    # Create weighted graph and cluster with Girvan-Newman algorithm, removing edges from lightest to heaviest
    graph = nx.from_pandas_edgelist(elusive_edges.to_pandas(), source="sample1", target="sample2", edge_attr=True)
    node_attributes = {k: {"read_size": v} for k,v in zip(read_size["sample"], read_size["read_size"])}
    nx.set_node_attributes(graph, node_attributes)

    def lightest(graph):
        u, v, _ = min(graph.edges(data="weight"), key=itemgetter(2))
        return (u, v)

    comp = community.girvan_newman(graph, most_valuable_edge=lightest)
    first_community = [tuple([c for c in components.connected_components(graph)])]

    communities = pl.DataFrame({"community": itertools.chain(first_community, comp)}
    ).lazy(
    ).with_columns(
        pl.col("community").apply(lambda x: [[i for i in s] for s in x])
    ).explode("community"
    ).with_columns(
        pl.col("community").arr.lengths().alias("length")
    ).with_columns(
        ((pl.col("length") <= MAX_RECOVERY_SAMPLES) &
        (pl.col("length") >= MIN_COASSEMBLY_SAMPLES)).alias("umbrella"),
        ((pl.col("length") >= MIN_COASSEMBLY_SAMPLES) &
        (pl.col("length") <= MAX_COASSEMBLY_SAMPLES)).alias("coassembly"),
        pl.col("community").apply(lambda x: [hash(y) for y in x]).hash().alias("coassembly_hash")
    ).filter(
        pl.col("coassembly") | pl.col("umbrella")
    ).unique(subset="coassembly_hash"
    ).collect()

    umbrellas = communities.filter(
        pl.col("umbrella")
    ).select(
        "community",
        pl.col("length").alias("umbrella_length"),
        pl.col("community").arr.sort().arr.join(",").alias("recover_samples")
    ).sort(
        "umbrella_length", descending=True
    )

    clusters = communities.filter(
        pl.col("coassembly")
    ).with_columns(
        pl.col("community").apply(
            lambda x: nx.subgraph_view(graph, filter_node = lambda y: y in x)
            ).alias("subgraphview"),
    ).with_columns(
        pl.col("subgraphview").apply(
            lambda x: sum([data["read_size"] for _,data in x.nodes(data=True)])
            ).alias("total_size"),
    ).filter(
        pl.when(MAX_COASSEMBLY_SIZE is not None)
        .then(pl.col("total_size") <= MAX_COASSEMBLY_SIZE)
        .otherwise(True)
    )

    if clusters.height == 0:
        return pl.DataFrame(schema=output_columns)

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
            (pl.col("recover_samples_list").arr.lengths().first() >= MAX_RECOVERY_SAMPLES) | (pl.col("relevant_edges").apply(lambda x: x.height) == 0)
        ).then(
            pl.Series("other_samples", [[]])
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

    elusive_edges = pl.read_csv(elusive_edges_path, sep="\t", dtypes={"target_ids": str})
    read_size = pl.read_csv(read_size_path, has_header=False, new_columns=["sample", "read_size"])

    clusters = pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES
        )
    clusters.write_csv(elusive_clusters_path, sep="\t")
