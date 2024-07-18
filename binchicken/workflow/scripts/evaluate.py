###################
### evaluate.py ###
###################
# Author: Samuel Aroney

import polars as pl
import os

SINGLEM_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
}

TARGET_COLUMNS = SINGLEM_COLUMNS | {
    "target": int,
}
APPRAISE_COLUMNS = SINGLEM_COLUMNS | {
    "found_in": str,
}

CLUSTER_COLUMNS = {
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": float,
    "recover_samples": str,
    "coassembly": str,
}
EDGE_COLUMNS = {
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
}

OUTPUT_COLUMNS = {
    "coassembly": str,
    "gene": str,
    "sequence": str,
    "genome": str,
    "target": str,
    "found_in": str,
    "source_samples": str,
    "source_num_hits": int,
    "source_coverage": float,
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

def evaluate(target_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table, recovered_bins):

    print(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(recovered_otu_table) == 0:
        empty_output = pl.DataFrame(schema=OUTPUT_COLUMNS)
        return empty_output, empty_output, pl.DataFrame(schema=SUMMARY_COLUMNS)

    # Load elusive clusters and edges (to match targets to coassemblies, with duplicates)
    sample_coassemblies = (
        elusive_clusters
        .select(pl.col("samples").str.split(","), "coassembly")
        .explode("samples")
    )

    elusive_edges = (
        elusive_edges
        .with_columns(
            pl.col("samples").str.split(",").list.eval(
                        pl.when(pl.all_horizontal(pl.element().is_in(sample_coassemblies.get_column("samples"))))
                        .then(pl.element())
                        .otherwise(pl.element().str.replace(r"(_|\.)1$", ""))
                        ),
            pl.col("samples").hash().alias("samples_hash")
            )
        .explode("samples")
        .join(sample_coassemblies, on="samples")
        .group_by("samples_hash", "coassembly")
        .agg(
            pl.first("target_ids"),
            pl.len().alias("count"),
            pl.first("cluster_size"),
            pl.col("samples").unique(),
            )
        .filter(pl.col("count") == pl.col("cluster_size"))
        .with_columns(target = pl.col("target_ids").str.split(","))
        .explode("target")
    )

    # Create otu table with original sequence, samples present, cluster id, target id and associated coassemblies
    sample_edges = (
        elusive_edges
        .group_by("coassembly", "target")
        .agg(
            source_samples = pl.col("samples")
                .flatten()
                .unique()
            )
        .explode("source_samples")
    )

    elusive_otu_table = (
        target_otu_table
        .with_columns(
            pl.col("target").cast(str),
            pl.lit(None).cast(str).alias("found_in"),
            source_samples =
                pl.when(pl.col("sample").is_in(sample_coassemblies.get_column("samples")))
                    .then(pl.col("sample"))
                .otherwise(pl.col("sample").str.replace(r"(_|\.)1$", "")),
            )
        .join(sample_edges, how="inner", on=["target", "source_samples"])
        .group_by("gene", "sequence", "coassembly", "taxonomy", "found_in", "target")
        .agg(
            pl.col("source_samples").sort().str.concat(","),
            source_num_hits = pl.sum("num_hits"),
            source_coverage = pl.sum("coverage"),
            )
    )

    # Add binned otu table to above with target NA
    nontarget_otu_table = (
        pl.concat([
            binned_otu_table
                .with_columns(target = pl.lit(None).cast(str)),
            target_otu_table
                .join(elusive_otu_table, on=["gene", "sequence"], how="anti")
                .drop("target")
                .with_columns(
                    found_in = pl.lit(None).cast(str),
                    target = pl.lit("None"),
                    ),
            ])
        .select([
            pl.col("sample").str.replace(r"_1$", "").str.replace(r"\.1$", ""),
            "gene", "sequence", "taxonomy", "found_in", "num_hits", "coverage", "target"
            ])
        .join(sample_coassemblies, left_on="sample", right_on="samples", how="left", coalesce=True)
        .drop_nulls("coassembly")
        .group_by(["gene", "sequence", "coassembly"])
        .agg([
            pl.first("taxonomy"),
            pl.first("found_in"),
            pl.first("target"),
            pl.col("sample").unique().sort().str.concat(",").alias("source_samples"),
            pl.sum("num_hits").alias("source_num_hits"),
            pl.sum("coverage").alias("source_coverage"),
            ])
        .unique()
    )

    haystack_otu_table = pl.concat([elusive_otu_table, nontarget_otu_table])

    # Load recovered otu table and join elusive and nontarget sequences
    recovered_otu_table = (
        recovered_otu_table
        .with_columns(
            pl.col("sample")
                .str.split("-")
                .list.to_struct(upper_bound=2, fields=["coassembly", "genome"])
            )
        .unnest("sample")
    )

    combined_otu_table = (
        recovered_otu_table
        .join(
            haystack_otu_table, on=["coassembly", "gene", "sequence"], how="outer_coalesce", suffix="old"
            )
        .select(
            "coassembly", "gene", "sequence", "genome", "target", "found_in",
            "source_samples", "source_num_hits", "source_coverage",
            pl.when(pl.col("taxonomy").is_null())
                .then(pl.col("taxonomyold"))
                .otherwise(pl.col("taxonomy"))
                .alias("taxonomy"),
            )
        .filter(
            (pl.col("coassembly").is_in(pl.lit(recovered_otu_table["coassembly"]))) &
            # Choose sequences where genome is present (from recovered) and/or target is present (from unbinned targets)
            ((pl.col("genome").is_not_null()) | (pl.col("target").is_not_null()))
            )
    )

    matches = (
        combined_otu_table
        .filter(
            ~pl.all_horizontal(pl.col(["target", "found_in", "source_samples"]).is_null())
            )
        .with_columns(
            target = pl.when(pl.col("target") == "None").then(pl.lit(None)).otherwise(pl.col("target"))
            )
    )

    unmatched = (
        combined_otu_table
        .filter(
            pl.all_horizontal(pl.col(["target", "found_in", "source_samples"]).is_null())
            )
    )

    # Summarise recovery stats
    summary = summarise_stats(matches, combined_otu_table, recovered_bins)

    return matches, unmatched, summary

def summarise_stats(matches, combined_otu_table, recovered_bins):
    recovered_hits = (
        combined_otu_table
        .filter(
            pl.col("genome").is_not_null()
            )
        .with_columns(
            target = pl.when(pl.col("target") == "None").then(pl.lit(None)).otherwise(pl.col("target"))
            )
    )

    summary = (
        matches
        .filter(pl.col("target").is_not_null())
        .with_columns(
            pl.when(pl.col("genome").is_not_null())
                .then(pl.lit("match"))
                .otherwise(pl.lit("nonmatch"))
                .alias("status")
            )
        .group_by([
            "coassembly", "status"
            ])
        .agg(
            pl.col("sequence").unique().len().alias("sequences"),
            pl.col("genome").unique().is_not_null().sum().alias("bins"),
            pl.col("taxonomy").unique().is_not_null().sum().alias("taxonomy"),
            )
        .join(
            # Duplicate sequences are counted multiple times to give a proportion at bin level
            recovered_hits
                .with_columns(
                    pl.when(
                        pl.col("found_in").is_not_null()
                    ).then(pl.lit("match")).otherwise(pl.lit("nonmatch")).alias("status")
                    )
                .group_by(["coassembly", "status"])
                .agg(pl.col("sequence").len().alias("nontarget_bin_sequences")),
            on=["coassembly", "status"], how="outer_coalesce"
            )
        .join(
            # Duplicate sequences are counted multiple times to give a proportion at bin level
            recovered_hits
                .with_columns(
                    pl.when(
                        pl.all_horizontal(pl.col(["target", "found_in"]).is_null()) & (pl.col("source_samples").is_not_null())
                    ).then(pl.lit("match")).otherwise(pl.lit("nonmatch")).alias("status")
                    )
                .group_by(["coassembly", "status"])
                .agg(pl.col("sequence").len().alias("nontarget_unbin_sequences")),
            on=["coassembly", "status"], how="outer_coalesce"
            )
        .join(
            # Duplicate sequences are counted multiple times to give a proportion at bin level
            recovered_hits
                .with_columns(
                    pl.when(
                        pl.all_horizontal(pl.col(["target", "found_in", "source_samples"]).is_null())
                    ).then(pl.lit("match")).otherwise(pl.lit("nonmatch")).alias("status")
                    )
                .group_by(["coassembly", "status"])
                .agg(pl.col("sequence").len().alias("novel_sequences")),
            on=["coassembly", "status"], how="outer_coalesce"
            )
        .melt(id_vars=["coassembly", "status"], variable_name="statistic", value_name="value")
        .fill_null(0)
    )

    recovered_counts = (
        pl.from_dict(recovered_bins)
        .melt(variable_name="name")
        .with_columns(
            pl.col("name").str.split("-").list.to_struct(upper_bound=2, fields=["coassembly", "genome"])
            )
        .unnest("name")
        .group_by("coassembly")
        .agg(pl.len().alias("n"))
        .with_columns(
            pl.lit("match").alias("status"),
            pl.lit("bins").alias("statistic"),
            )
        .join(summary.rename({"value": "match"}), on=["coassembly", "status", "statistic"], how="left", coalesce=True)
        .select(
            "coassembly",
            pl.lit("nonmatch").alias("status"),
            "statistic",
            pl.col("n").alias("value") - pl.col("match")
            )
    )

    summary = (
        pl.concat([
            summary.filter((pl.col("statistic") != "bins") | (pl.col("status") != "nonmatch")),
            recovered_counts
            ])
        .pivot(index=["coassembly", "statistic"], values="value", columns="status", aggregate_function=None)
        .select(
            "coassembly", "statistic",
            pl.when(pl.col("statistic").is_in(["sequences", "taxonomy"])).then(pl.lit("targets")).otherwise(pl.lit("recovery")).alias("within"),
            "match", "nonmatch",
            pl.col("match").alias("total") + pl.col("nonmatch"),
            (pl.col("match") / (pl.col("match") + pl.col("nonmatch")) * 100).round(2).alias("match_percent"),
            )
        .sort(["coassembly", "within", "statistic"], descending=[False, True, False])
    )

    return summary

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    target_path = snakemake.params.target_otu_table
    binned_path = snakemake.params.binned_otu_table
    elusive_clusters_path = snakemake.params.elusive_clusters
    elusive_edges_path = snakemake.params.elusive_edges
    recovered_otu_table_path = snakemake.input.recovered_otu_table
    recovered_bins = snakemake.params.recovered_bins
    matched_hits_path = snakemake.output.matched_hits
    novel_hits_path = snakemake.output.novel_hits
    summary_stats_path = snakemake.output.summary_stats

    target_otu_table = pl.read_csv(target_path, separator="\t", dtypes=TARGET_COLUMNS)
    binned_otu_table = pl.read_csv(binned_path, separator="\t", dtypes=APPRAISE_COLUMNS)
    elusive_clusters = pl.read_csv(elusive_clusters_path, separator="\t", dtypes=CLUSTER_COLUMNS)
    elusive_edges = pl.read_csv(elusive_edges_path, separator="\t", dtypes=EDGE_COLUMNS)
    recovered_otu_table = pl.read_csv(recovered_otu_table_path, separator="\t", dtypes=SINGLEM_COLUMNS)

    matches, unmatched, summary = evaluate(target_otu_table, binned_otu_table, elusive_clusters, elusive_edges, recovered_otu_table, recovered_bins)
    # Export hits matching elusive targets
    matches.write_csv(matched_hits_path, separator="\t")
    # Export non-elusive sequence hits
    unmatched.write_csv(novel_hits_path, separator="\t")
    # Export summary stats
    summary.write_csv(summary_stats_path, separator="\t")
