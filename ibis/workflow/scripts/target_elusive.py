#########################
### target_elusive.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os
import logging

EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
    }

def pipeline(
    unbinned,
    MIN_COASSEMBLY_COVERAGE=10,
    TAXA_OF_INTEREST="",
    MAX_COASSEMBLY_SAMPLES=2):

    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    if len(unbinned) == 0:
        logging.warning("No unbinned sequences found")
        return unbinned.rename({"found_in": "target"}), pl.DataFrame(schema=EDGES_COLUMNS)

    if MAX_COASSEMBLY_SAMPLES < 2:
        # Set to 2 to produce paired edges
        MAX_COASSEMBLY_SAMPLES = 2

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    logging.info("Grouping hits by marker gene sequences to form targets")
    unbinned = (
        unbinned
        .drop("found_in")
        .with_row_count("target")
        .select(
            "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
            pl.first("target").over(["gene", "sequence"]).rank("dense") - 1,
            )
        .with_columns(
            pl.col("target").cast(pl.Utf8)
            )
    )

    def process_groups(df):
        if df.height == 1:
            return pl.DataFrame(schema={"style": str, "cluster_size": pl.Int64, "samples": str, "target": str})

        # Direct matching samples in pairs with coverage > MIN_COASSEMBLY_COVERAGE
        dfs = [
            df
            .join(df, how="cross", suffix="_2")
            .filter(
                pl.col("samples").list.last().str.encode("hex") < pl.col("samples_2").list.first().str.encode("hex")
                )
            .select([
                "target",
                pl.col("coverage") + pl.col("coverage_2").alias("coverage"),
                pl.col("samples").list.concat(pl.col("samples_2")),
                ])
            .filter(pl.col("samples").list.lengths() > 1)
            .filter(pl.col("coverage") > MIN_COASSEMBLY_COVERAGE)
            .select(
                pl.lit("match").alias("style"),
                pl.lit(2).cast(pl.Int64).alias("cluster_size"),
                pl.col("samples").list.join(","),
                pl.col("target"),
                )
        ]

        # Pool samples with coverage > MIN_COASSEMBLY_COVERAGE / N for clusters of size N: 3 to MAX_COASSEMBLY_SAMPLES
        cluster_sizes = pl.DataFrame({"cluster_size": range(3, MAX_COASSEMBLY_SAMPLES+1)})
        dfs.append(
            df
            .join(cluster_sizes, how="cross")
            .filter(pl.col("coverage") > float(MIN_COASSEMBLY_COVERAGE) / pl.col("cluster_size").cast(float))
            .groupby("target", "cluster_size")
            .agg(pl.concat_list("samples").flatten())
            .filter(pl.col("samples").list.lengths() >= pl.col("cluster_size"))
            .select(
                pl.lit("pool").alias("style"),
                pl.col("cluster_size"),
                pl.col("samples").list.sort().list.join(","),
                pl.col("target"),
                )
        )

        return pl.concat(dfs)

    logging.info("Grouping targets into paired matches and pooled samples for clusters of size 3+")
    sparse_edges = (
        unbinned
        .select(
            "target",
            "coverage",
            samples = pl.col("sample").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)),
            )
        .groupby("target")
        .apply(process_groups)
        .groupby(["style", "cluster_size", "samples"])
        .agg(target_ids = pl.col("target").sort().str.concat(","))
    )

    return unbinned, sparse_edges

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    MIN_COASSEMBLY_COVERAGE = snakemake.params.min_coassembly_coverage
    MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    unbinned_path = snakemake.input.unbinned
    targets_path = snakemake.output.output_targets
    edges_path = snakemake.output.output_edges

    unbinned = pl.read_csv(unbinned_path, separator="\t")

    targets, edges = pipeline(
        unbinned,
        MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
        TAXA_OF_INTEREST=TAXA_OF_INTEREST,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        )
    targets.write_csv(targets_path, separator="\t")
    edges.write_csv(edges_path, separator="\t")
