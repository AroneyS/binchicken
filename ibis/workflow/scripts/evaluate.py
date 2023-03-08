###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import polars as pl
import os

OUTPUT_COLUMNS={
    "coassembly": str,
    "gene": str,
    "sequence": str,
    "genome": str,
    "target": str,
    "found_in": str,
    "taxonomy": str,
    }
SUMMARY_COLUMNS = {
    "coassembly": str,
    "statistic": str,
    "within": str,
    "match": int,
    "nonmatch": int,
    "total": int,
    "match_percent": float,
    }

def evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table, recovered_bins):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    if len(recovered_otu_table) == 0:
        empty_output = pl.DataFrame(schema=OUTPUT_COLUMNS)
        return empty_output, empty_output, pl.DataFrame(schema=SUMMARY_COLUMNS)

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

    # Summarise recovery stats
    summary = summarise_stats(matches, combined_otu_table, recovered_bins)

    return matches, unmatched, summary

def summarise_stats(matches, combined_otu_table, recovered_bins):
    recovered_hits = combined_otu_table.filter(
        pl.col("genome").is_not_null()
    )

    summary = matches.filter(
        pl.col("target").is_not_null()
    ).with_columns(
        pl.when(pl.col("genome").is_not_null())
        .then("match")
        .otherwise("nonmatch")
        .alias("status")
    ).groupby([
        "coassembly", "status"
    ]).agg(
        pl.col("sequence").unique().len().alias("sequences"),
        pl.col("genome").unique().is_not_null().sum().alias("bins"),
        pl.col("taxonomy").unique().is_not_null().sum().alias("taxonomy"),
    ).join(
        # Duplicate sequences are counted multiple times to give a proportion at bin level
        recovered_hits.with_columns(
            pl.when(
                pl.col("found_in").is_not_null()
            ).then("match").otherwise("nonmatch").alias("status")
        ).groupby([
            "coassembly", "status"
        ]).agg(
            pl.col("sequence").len().alias("nontarget_sequences")
        ),
        on=["coassembly", "status"], how="outer"
    ).join(
        # Duplicate sequences are counted multiple times to give a proportion at bin level
        recovered_hits.with_columns(
            pl.when(
                (pl.col("found_in").is_null()) & (pl.col("target").is_null())
            ).then("match").otherwise("nonmatch").alias("status")
        ).groupby([
            "coassembly", "status"
        ]).agg(
            pl.col("sequence").len().alias("novel_sequences")
        ),
        on=["coassembly", "status"], how="outer"
    ).melt(
        id_vars=["coassembly", "status"], variable_name="statistic", value_name="value"
    ).fill_null(
        0
    )

    recovered_counts = pl.from_dict(recovered_bins
    ).melt(variable_name="name"
    ).with_columns(
        pl.col("name").str.split("-").arr.to_struct(upper_bound=2, name_generator=lambda id: ["coassembly", "genome"][id])
    ).unnest("name"
    ).groupby("coassembly"
    ).agg(
        pl.count().alias("n")
    ).with_columns(
        pl.lit("match").alias("status"),
        pl.lit("bins").alias("statistic"),
    ).join(summary.rename({"value": "match"}), on=["coassembly", "status", "statistic"], how="left"
    ).select(
        "coassembly",
        pl.lit("nonmatch").alias("status"),
        "statistic",
        pl.col("n").alias("value") - pl.col("match")
    )

    summary = pl.concat([
        summary.filter((pl.col("statistic") != "bins") | (pl.col("status") != "nonmatch")),
        recovered_counts
    ]).pivot(
        index=["coassembly", "statistic"], values="value", columns="status"
    ).select(
        "coassembly", "statistic",
        pl.when(pl.col("statistic").is_in(["sequences", "taxonomy"])).then("targets").otherwise("recovery").alias("within"),
        "match", "nonmatch",
        pl.col("match").alias("total") + pl.col("nonmatch"),
        (pl.col("match") / (pl.col("match") + pl.col("nonmatch")) * 100).round(2).alias("match_percent"),
    ).sort(
        ["coassembly", "within", "statistic"], descending=[False, True, False]
    )

    return summary

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    unbinned_path = snakemake.params.unbinned_otu_table
    binned_path = snakemake.params.binned_otu_table
    elusive_clusters_path = snakemake.params.elusive_clusters
    elusive_edges_path = snakemake.params.elusive_edges
    recovered_otu_table_path = snakemake.input.recovered_otu_table
    recovered_bins = snakemake.params.recovered_bins
    matched_hits_path = snakemake.output.matched_hits
    novel_hits_path = snakemake.output.novel_hits
    summary_stats_path = snakemake.output.summary_stats

    unbinned_otu_table = pl.read_csv(unbinned_path, sep="\t")
    binned_otu_table = pl.read_csv(binned_path, sep="\t")
    elusive_clusters = pl.read_csv(elusive_clusters_path, sep="\t")
    elusive_edges = pl.read_csv(elusive_edges_path, sep="\t")
    recovered_otu_table = pl.read_csv(recovered_otu_table_path, sep="\t")

    matches, unmatched, summary = evaluate(unbinned_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table, recovered_bins)
    # Export hits matching elusive targets
    matches.write_csv(matched_hits_path, sep="\t")
    # Export non-elusive sequence hits
    unmatched.write_csv(novel_hits_path, sep="\t")
    # Export summary stats
    summary.write_csv(summary_stats_path, sep="\t")
