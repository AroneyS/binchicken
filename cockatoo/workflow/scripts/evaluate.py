###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import pandas as pd
import re

def evaluate(unbinned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table):
    # Load otu table of unbinned sequences and get unique id for each sequence (to match sequences to target id)
    unbinned_otu_table = unbinned_otu_table.groupby(["gene", "sequence"]).first()[["target"]].reset_index()
    unbinned_otu_table.dropna(inplace=True)
    unbinned_otu_table["target"] = unbinned_otu_table["target"].astype(str)

    # Load elusive clusters and edges (to match targets to coassemblies, with duplicates)
    elusive_clusters["samples"] = elusive_clusters["samples"].str.split(",")
    elusive_clusters = elusive_clusters.explode("samples").set_index("samples")[["coassembly"]]

    elusive_edges["sample1"] = elusive_edges["sample1"].apply(lambda x: re.sub("\.1$", "", x))
    elusive_edges["sample2"] = elusive_edges["sample2"].apply(lambda x: re.sub("\.1$", "", x))
    elusive_edges = (elusive_edges
        .set_index("sample1")
        .join(elusive_clusters)
        .reset_index()
        .rename(columns={"index": "sample1"})
        .set_index("sample2")
        .join(elusive_clusters, rsuffix="2")
        .reset_index()
        .rename(columns={"index": "sample2"})
        )

    coassembly_edges = elusive_edges[elusive_edges["coassembly"] == elusive_edges["coassembly2"]].copy()
    coassembly_edges["target"] = coassembly_edges["target_ids"].str.split(",")
    coassembly_edges = coassembly_edges.explode("target").drop_duplicates(["target", "coassembly"]).set_index("target")[["coassembly"]]

    # Create otu table with original sequence, cluster id, target id and associated coassemblies
    elusive_otu_table = coassembly_edges.join(unbinned_otu_table.set_index("target")).reset_index().set_index(["coassembly", "gene", "sequence"])

    # Load recovered otu table and join elusive sequences
    recovered_otu_table[["coassembly", "genome"]] = recovered_otu_table["sample"].str.split("-", 1, expand=True)
    recovered_coassemblies = set(recovered_otu_table["coassembly"])

    combined_otu_table = recovered_otu_table.set_index(["coassembly", "gene", "sequence"])[["genome", "taxonomy"]].join(elusive_otu_table, how="outer").reset_index()
    combined_otu_table.insert(len(combined_otu_table.columns)-1, "taxonomy", combined_otu_table.pop("taxonomy"))
    combined_otu_table = combined_otu_table[combined_otu_table["coassembly"].isin(recovered_coassemblies)]
    matches = combined_otu_table.dropna(subset = ["target"]).reset_index(drop=True)
    unmatched = combined_otu_table[combined_otu_table["target"].isnull()].reset_index(drop=True)

    return matches, unmatched


if __name__ == "__main__":
    unbinned_path = snakemake.params.unbinned_otu_table
    elusive_clusters_path = snakemake.params.elusive_clusters
    elusive_edges_path = snakemake.params.elusive_edges
    recovered_otu_table_path = snakemake.input.recovered_otu_table
    unbinned_hits_path = snakemake.output.unbinned_hits
    novel_hits_path = snakemake.output.novel_hits

    unbinned_otu_table = pd.read_csv(unbinned_path, sep="\t")
    elusive_clusters = pd.read_csv(elusive_clusters_path, sep="\t")
    elusive_edges = pd.read_csv(elusive_edges_path, sep="\t")
    recovered_otu_table = pd.read_csv(recovered_otu_table_path, sep="\t")

    matches, unmatched = evaluate(unbinned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table)
    # Export hits matching elusive targets
    matches.to_csv(unbinned_hits_path, sep="\t", index=False)
    # Export non-elusive sequence hits
    unmatched.to_csv(novel_hits_path, sep="\t", index=False)
