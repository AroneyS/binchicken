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
    MIN_COASSEMBLY_COVERAGE = 10,
    TAXA_OF_INTEREST = ""):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    if len(unbinned) == 0:
        return unbinned.rename({"found_in": "target"}), pl.DataFrame(schema=EDGES_COLUMNS)

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    # Group hits by sequence within genes and number to form targets
    unbinned = unbinned.drop("found_in"
    ).with_row_count("target"
    ).select(
        "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
        pl.first("target").over(["gene", "sequence"]).rank("dense") - 1,
    ).with_columns(
        pl.col("target").cast(pl.Utf8)
    )

    # Within each target cluster, find pairs of samples with combined coverage > MIN_COVERAGE
    sample_pairs = unbinned.lazy().join(unbinned.lazy(), on="target", how="left", suffix="_2"
    ).filter(
        (pl.col("sample").str.encode("hex") < pl.col("sample_2").str.encode("hex")) &
        (pl.col("coverage") + pl.col("coverage_2") > MIN_COASSEMBLY_COVERAGE)
    ).collect(streaming=True)

    # Create weighted graph with nodes as samples and edges weighted by the number of clusters supporting that co-assembly (networkx)
    # Output sparse matrix with sample1/sample2/edge weight (number of supporting clusters with coverage > threshold)
    sparse_edges = sample_pairs.groupby([
        "sample", "sample_2"
    ]).agg([
        pl.count().alias("weight").cast(int),
        pl.col("target").sort().str.concat(",").alias("target_ids")
    ]).select([
        pl.concat_str(["sample", "sample_2"], separator=",").alias("samples"),
        "weight", "target_ids",
    ])

    return unbinned, sparse_edges

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    MIN_COASSEMBLY_COVERAGE = snakemake.params.min_coassembly_coverage
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    unbinned_path = snakemake.input.unbinned
    targets_path = snakemake.output.output_targets
    edges_path = snakemake.output.output_edges

    unbinned = pl.read_csv(unbinned_path, separator="\t")

    targets, edges = pipeline(
        unbinned,
        MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
        TAXA_OF_INTEREST=TAXA_OF_INTEREST
        )
    targets.write_csv(targets_path, separator="\t")
    edges.write_csv(edges_path, separator="\t")
