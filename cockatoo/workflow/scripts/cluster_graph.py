########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import pandas as pd
import re
import networkx as nx
from networkx.algorithms import community
from operator import itemgetter

MAX_COASSEMBLY_SIZE = snakemake.params.max_coassembly_size * 10**9 if snakemake.params.max_coassembly_size else None
MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
MIN_COASSEMBLY_SAMPLES = snakemake.params.num_coassembly_samples
MAX_RECOVERY_SAMPLES = snakemake.params.max_recovery_samples

################
### Pipeline ###
################
# Load data, fixing sample names and preparing read_size
elusive_edges = pd.read_csv(snakemake.input.elusive_edges, sep="\t")
elusive_edges["sample1"] = elusive_edges["sample1"].apply(lambda x: re.sub("\.1$", "", x))
elusive_edges["sample2"] = elusive_edges["sample2"].apply(lambda x: re.sub("\.1$", "", x))

read_size = pd.read_csv(snakemake.input.read_size, names = ["sample", "read_size"])
node_attributes = {k: {"read_size": v} for k,v in zip(read_size["sample"], read_size["read_size"])}

# Create weighted graph and cluster with Girvan-Newman algorithm, removing edges from lightest to heaviest
graph = nx.from_pandas_edgelist(elusive_edges, source="sample1", target="sample2", edge_attr=True)
nx.set_node_attributes(graph, node_attributes)

def lightest(graph):
    u, v, _ = min(graph.edges(data="weight"), key=itemgetter(2))
    return (u, v)

comp = community.girvan_newman(graph, most_valuable_edge=lightest)
clusters = pd.DataFrame(columns=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples"])
umbrellas = []
for iteration in comp:
    communities = tuple(sorted(c) for c in iteration)
    for c in communities:
        if len(c) <= MAX_RECOVERY_SAMPLES and len(c) >= MIN_COASSEMBLY_SAMPLES and not c in umbrellas:
            umbrellas.append(c)
        if len(c) > MAX_COASSEMBLY_SAMPLES: continue
        if len(c) < MIN_COASSEMBLY_SAMPLES: continue

        subgraphview = nx.subgraph_view(
            graph,
            filter_node = lambda x: x in c
            )

        total_size = sum([data["read_size"] for _,data in subgraphview.nodes(data=True)])
        if MAX_COASSEMBLY_SIZE and total_size > MAX_COASSEMBLY_SIZE: continue

        edgeview = subgraphview.edges(data=True)
        total_weight = sum([data["weight"] for _,_,data in edgeview])
        total_targets = len(set([i for l in [data["target_ids"].split(",") for _,_,data in edgeview] for i in l]))

        suitable_indices = [i for i,u in enumerate(umbrellas) if all([s in u for s in c])]
        suitable_umbrellas = [umbrellas[i] for i in suitable_indices]
        best_umbrella = suitable_umbrellas[0]

        df = pd.DataFrame({
            "samples": ",".join(c),
            "length": len(c),
            "total_weight": total_weight,
            "total_targets": total_targets,
            "total_size": total_size,
            "recover_samples": ",".join(best_umbrella)
            },
            index = [0])
        clusters = pd.concat([clusters, df], ignore_index=True)

clusters.drop_duplicates(inplace=True)
if MAX_COASSEMBLY_SAMPLES == 1:
    def find_top_samples(sample):
        sample_edges = elusive_edges[(elusive_edges["sample1"] == sample) | (elusive_edges["sample2"] == sample)].copy()
        sample_edges.sort_values(by="weight", ascending=False, inplace=True)
        recover_samples = sample_edges.head(MAX_RECOVERY_SAMPLES - 1)
        recover_samples = recover_samples["sample1"].to_list() + recover_samples["sample2"].to_list()
        return ",".join(sorted(set(recover_samples)))

    clusters["recover_samples"] = clusters["samples"].apply(find_top_samples)

clusters.sort_values(by=["total_targets", "samples"], ascending=False, inplace=True)
clusters.reset_index(drop=True, inplace=True)
clusters["coassembly"] = (clusters.reset_index()["index"].apply(lambda x: "coassembly_" + str(x)))
clusters.to_csv(snakemake.output.elusive_clusters, sep="\t", index=False)
