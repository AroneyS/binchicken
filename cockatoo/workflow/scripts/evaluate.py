###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import pandas as pd
import re

# Load otu table of unbinned sequences and get unique id for each sequence (to match sequences to target id)
unbinned_otu_table = pd.read_csv(snakemake.params.unbinned_otu_table, sep="\t")
unbinned_otu_table = unbinned_otu_table.groupby(["gene", "sequence"]).first()[["target", "taxonomy"]].reset_index()
unbinned_otu_table.dropna(inplace=True)
unbinned_otu_table["target"] = unbinned_otu_table["target"].astype(str)

# Load elusive clusters and edges (to match targets to coassemblies, with duplicates)
elusive_clusters = pd.read_csv(snakemake.params.elusive_clusters, sep="\t")
elusive_clusters["samples"] = elusive_clusters["samples"].str.split(",")
elusive_clusters = elusive_clusters.explode("samples").set_index("samples")[["coassembly"]]

elusive_edges = pd.read_csv(snakemake.params.elusive_edges, sep="\t")
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
recovered_otu_table = pd.read_csv(snakemake.input.recovered_otu_table, sep="\t")
recovered_otu_table[["coassembly", "genome"]] = recovered_otu_table["sample"].str.split("-", 1, expand=True)
recovered_coassemblies = set(recovered_otu_table["coassembly"])

combined_otu_table = recovered_otu_table.set_index(["coassembly", "gene", "sequence"])[["genome"]].join(elusive_otu_table, how="outer").reset_index()
combined_otu_table = combined_otu_table[combined_otu_table["coassembly"].isin(recovered_coassemblies)]
matches = combined_otu_table.dropna(subset = ["target"])
unmatched = combined_otu_table[combined_otu_table["target"].isnull()]

# Export hits matching elusive targets
matches.to_csv(snakemake.output.unbinned_hits, sep="\t", index=False)
# Export non-elusive sequence hits
unmatched.to_csv(snakemake.output.novel_hits, sep="\t", index=False)
