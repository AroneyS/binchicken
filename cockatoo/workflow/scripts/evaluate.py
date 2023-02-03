###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import pandas as pd
import re

def evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table):

    output_columns = ["coassembly", "gene", "sequence", "genome", "target", "found_in", "taxonomy"]

    if len(recovered_otu_table) == 0:
        empty_output = pd.DataFrame(columns=output_columns)
        return empty_output, empty_output

    # Load otu table of unbinned sequences and get unique id for each sequence (to match sequences to target id)
    unbinned_otu_table = unbinned_otu_table.groupby(["gene", "sequence"]).first()[["target", "taxonomy"]].reset_index()
    unbinned_otu_table.dropna(inplace=True)
    unbinned_otu_table["target"] = unbinned_otu_table["target"].astype(str)

    # Load elusive clusters and edges (to match targets to coassemblies, with duplicates)
    elusive_clusters["samples"] = elusive_clusters["samples"].str.split(",")
    elusive_clusters = elusive_clusters.explode("samples").set_index("samples")[["coassembly"]]

    elusive_edges["sample1"] = elusive_edges["sample1"].apply(lambda x: re.sub(r"\.1$", "", x))
    elusive_edges["sample2"] = elusive_edges["sample2"].apply(lambda x: re.sub(r"\.1$", "", x))
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
    elusive_otu_table["found_in"] = None

    # Add binned otu table to above with target NA
    binned_otu_table["sample"] = binned_otu_table["sample"].apply(lambda x: re.sub(r"\.1$", "", x))
    nontarget_otu_table = (binned_otu_table
        .set_index(["sample"])[["gene", "sequence", "taxonomy", "found_in"]]
        .join(elusive_clusters)
        .dropna(subset=["coassembly"])
        .set_index(["coassembly", "gene", "sequence"])
        .drop_duplicates()
        )
    nontarget_otu_table["target"] = None

    haystack_otu_table = pd.concat([elusive_otu_table, nontarget_otu_table])

    # Load recovered otu table and join elusive and nontarget sequences
    recovered_otu_table[["coassembly", "genome"]] = recovered_otu_table["sample"].str.split("-", n=1, expand=True)
    recovered_coassemblies = set(recovered_otu_table["coassembly"])

    combined_otu_table = recovered_otu_table.set_index(["coassembly", "gene", "sequence"])[["genome", "taxonomy"]].join(haystack_otu_table, how="outer", rsuffix="old").reset_index()
    combined_otu_table["taxonomy"] = combined_otu_table["taxonomy"].combine(combined_otu_table["taxonomyold"], lambda a,b: a if not pd.isna(a) else b)
    combined_otu_table = combined_otu_table.drop("taxonomyold", axis=1)
    combined_otu_table.insert(len(combined_otu_table.columns)-1, "taxonomy", combined_otu_table.pop("taxonomy"))
    combined_otu_table = combined_otu_table[combined_otu_table["coassembly"].isin(recovered_coassemblies)]
    # Choose sequences where genome is present (from recovered) and/or target is present (from unbinned targets)
    combined_otu_table = combined_otu_table[combined_otu_table["genome"].notnull() | combined_otu_table["target"].notnull()]

    matches = combined_otu_table.dropna(subset=["target", "found_in"], how="all").reset_index(drop=True)
    unmatched = combined_otu_table[combined_otu_table["target"].isnull() & combined_otu_table["found_in"].isnull()].reset_index(drop=True)

    return matches, unmatched


if __name__ == "__main__":
    unbinned_path = snakemake.params.unbinned_otu_table
    binned_path = snakemake.params.binned_otu_table
    elusive_clusters_path = snakemake.params.elusive_clusters
    elusive_edges_path = snakemake.params.elusive_edges
    recovered_otu_table_path = snakemake.input.recovered_otu_table
    matched_hits_path = snakemake.output.matched_hits
    novel_hits_path = snakemake.output.novel_hits

    unbinned_otu_table = pd.read_csv(unbinned_path, sep="\t")
    binned_otu_table = pd.read_csv(binned_path, sep="\t")
    elusive_clusters = pd.read_csv(elusive_clusters_path, sep="\t")
    elusive_edges = pd.read_csv(elusive_edges_path, sep="\t")
    recovered_otu_table = pd.read_csv(recovered_otu_table_path, sep="\t")

    matches, unmatched = evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table)
    # Export hits matching elusive targets
    matches.to_csv(matched_hits_path, sep="\t", index=False)
    # Export non-elusive sequence hits
    unmatched.to_csv(novel_hits_path, sep="\t", index=False)
