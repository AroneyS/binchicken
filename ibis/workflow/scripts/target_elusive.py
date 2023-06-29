#########################
### target_elusive.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os

EDGES_COLUMNS={
    "samples": str,
    "weight": int,
    "target_ids": str,
    }

def pipeline(
    unbinned,
    MIN_COASSEMBLY_COVERAGE=10,
    TAXA_OF_INTEREST="",
    MAX_COASSEMBLY_SAMPLES=2):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    if len(unbinned) == 0:
        return unbinned.rename({"found_in": "target"}), pl.DataFrame(schema=EDGES_COLUMNS)

    if MAX_COASSEMBLY_SAMPLES < 2:
        # Set to 2 to produce paired edges
        MAX_COASSEMBLY_SAMPLES = 2

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    # Group hits by sequence within genes and number to form targets
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
        dfs = [df.lazy()]
        for _ in range(1, MAX_COASSEMBLY_SAMPLES):
            dfs.append(
                dfs[-1]
                .join(df.lazy(), on="target", how="left", suffix="_2")
                .filter(
                    (pl.col("samples").arr.contains(pl.col("samples_2").arr.first()).is_not()) &
                    (pl.col("samples").arr.last().str.encode("hex") < pl.col("samples_2").arr.first().str.encode("hex"))
                    )
                .select([
                    "target",
                    pl.col("coverage") + pl.col("coverage_2").alias("coverage"),
                    pl.col("samples").arr.concat(pl.col("samples_2")),
                    ])
            )

        filtered = (
            pl.concat(
                pl.collect_all(dfs)
                )
            .filter(pl.col("samples").arr.lengths() > 1)
            .filter(pl.col("coverage") > MIN_COASSEMBLY_COVERAGE)
            .with_columns(pl.col("samples").arr.join(","))
        )

        return filtered

    sparse_edges = (
        unbinned
        .select(
            "target",
            "coverage",
            samples = pl.col("sample").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)),
            )
        .groupby("target")
        .apply(process_groups)
        .groupby("samples", maintain_order=True)
        .agg(
            weight = pl.count().cast(int),
            target_ids = pl.col("target").sort().str.concat(","),
            )
    )

    return unbinned, sparse_edges

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

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
