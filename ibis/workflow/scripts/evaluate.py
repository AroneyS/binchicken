###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import polars as pl

OUTPUT_COLUMNS={
    "coassembly": str,
    "gene": str,
    "sequence": str,
    "genome": str,
    "target": str,
    "found_in": str,
    "taxonomy": str,
    }

def evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table):

    if len(recovered_otu_table) == 0:
        empty_output = pl.DataFrame(schema=OUTPUT_COLUMNS)
        return empty_output, empty_output

    # Load otu table of unbinned sequences and get unique id for each sequence (to match sequences to target id)
    unbinned_otu_table = unbinned_otu_table.select([
        "gene", "sequence",
        pl.first("target").over(["gene", "sequence"]).cast(str),
        pl.first("taxonomy").over(["gene", "sequence"]),
    ]).unique().drop_nulls()

    # Load elusive clusters and edges (to match targets to coassemblies, with duplicates)
    sample_coassemblies = elusive_clusters.select(
        pl.col("samples").str.split(","),
        "coassembly"
    ).explode("samples")

    coassembly_edges = elusive_edges.with_columns(
        pl.col("sample1").str.replace(r"\.1$", ""),
        pl.col("sample2").str.replace(r"\.1$", ""),
    ).join(
        sample_coassemblies, left_on="sample1", right_on="samples", how="left"
    ).join(
        sample_coassemblies, left_on="sample2", right_on="samples", how="left", suffix="2"
    ).filter(
        pl.col("coassembly") == pl.col("coassembly2")
    ).with_columns(
        pl.col("target_ids").str.split(",").alias("target")
    ).explode("target"
    ).select(["target", "coassembly"]
    ).unique()

    # Create otu table with original sequence, cluster id, target id and associated coassemblies
    elusive_otu_table = coassembly_edges.join(
        unbinned_otu_table, on="target", how="left"
    ).select(
        "gene", "sequence", "taxonomy",
        pl.lit(None).cast(str).alias("found_in"),
        "coassembly", "target",
    )

    # Add binned otu table to above with target NA
    nontarget_otu_table = binned_otu_table.select([
        pl.col("sample").str.replace(r"\.1$", ""),
        "gene", "sequence", "taxonomy", "found_in"
    ]).join(
        sample_coassemblies, left_on="sample", right_on="samples", how="left"
    ).drop("sample"
    ).drop_nulls("coassembly"
    ).unique(
    ).with_columns(
        pl.lit(None).cast(str).alias("target")
    )

    haystack_otu_table = pl.concat([elusive_otu_table, nontarget_otu_table])

    # Load recovered otu table and join elusive and nontarget sequences
    recovered_otu_table = recovered_otu_table.with_columns(
        pl.col("sample").str.split("-").arr.to_struct(upper_bound=2, name_generator=lambda id: ["coassembly", "genome"][id])
    ).unnest("sample"
    )

    combined_otu_table = recovered_otu_table.join(
        haystack_otu_table, on=["coassembly", "gene", "sequence"], how="outer", suffix="old"
    ).select(
        "coassembly", "gene", "sequence", "genome", "target", "found_in",
        pl.when(pl.col("taxonomy").is_null())
        .then(pl.col("taxonomyold"))
        .otherwise(pl.col("taxonomy"))
        .alias("taxonomy"),
    ).filter(
        (pl.col("coassembly").is_in(pl.lit(recovered_otu_table["coassembly"]))) &
        # Choose sequences where genome is present (from recovered) and/or target is present (from unbinned targets)
        ((pl.col("genome").is_not_null()) | (pl.col("target").is_not_null()))
    )

    matches = combined_otu_table.filter(
        ~pl.all(pl.col(["target", "found_in"]).is_null())
    )

    unmatched = combined_otu_table.filter(
        (pl.col("target").is_null()) & (pl.col("found_in").is_null())
    )

    return matches, unmatched


if __name__ == "__main__":
    unbinned_path = snakemake.params.unbinned_otu_table
    binned_path = snakemake.params.binned_otu_table
    elusive_clusters_path = snakemake.params.elusive_clusters
    elusive_edges_path = snakemake.params.elusive_edges
    recovered_otu_table_path = snakemake.input.recovered_otu_table
    matched_hits_path = snakemake.output.matched_hits
    novel_hits_path = snakemake.output.novel_hits

    unbinned_otu_table = pl.read_csv(unbinned_path, sep="\t")
    binned_otu_table = pl.read_csv(binned_path, sep="\t")
    elusive_clusters = pl.read_csv(elusive_clusters_path, sep="\t")
    elusive_edges = pl.read_csv(elusive_edges_path, sep="\t")
    recovered_otu_table = pl.read_csv(recovered_otu_table_path, sep="\t")

    matches, unmatched = evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table)
    # Export hits matching elusive targets
    matches.write_csv(matched_hits_path, sep="\t")
    # Export non-elusive sequence hits
    unmatched.write_csv(novel_hits_path, sep="\t")
