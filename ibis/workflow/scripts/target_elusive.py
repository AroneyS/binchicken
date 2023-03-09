#########################
### target_elusive.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os
from itertools import combinations

EDGES_COLUMNS={
    "taxa_group": str,
    "weight": int,
    "target_ids": str,
    "sample1": str,
    "sample2": str
    }

def pipeline(
    unbinned,
    MIN_COASSEMBLY_COVERAGE = 10,
    TAXA_OF_INTEREST = ""):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    if len(unbinned) == 0:
        return unbinned.rename({"found_in": "target"}), pl.DataFrame(schema=EDGES_COLUMNS)

    TAXA_LEVEL_SEP = "; "
    taxa_levels = {
        "d": 1,
        "p": 2,
        "c": 3,
        "o": 4,
        "f": 5,
        "g": 6,
        "s": 7,
    }
    try:
        TAXA_LEVEL_OF_INTEREST = taxa_levels[TAXA_OF_INTEREST[0]]
    except IndexError:
        TAXA_LEVEL_OF_INTEREST = 0

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
    ).with_columns(
        # Select group of interest (corresponding taxa level in brackets)
        # e.g. Root (0) or d__Bacteria (1) or p__Acidobacteriota (2) or c__Acidobacteriae (3) or o__Bryobacterales (4)
        # or f__Bryobacteraceae (5) or g__Solibacter (6) or s__Solibacter sp003136215 (7)
        pl.col("taxonomy").str.split(TAXA_LEVEL_SEP).arr.get(TAXA_LEVEL_OF_INTEREST).alias("taxa_group"),
    ).filter(
        (pl.col("taxa_group").is_not_null()) &
        ((TAXA_OF_INTEREST == "") |
        (pl.col("taxa_group") == TAXA_OF_INTEREST))
    ).collect(streaming=True)

    # Create weighted graph with nodes as samples and edges weighted by the number of clusters supporting that co-assembly (networkx)
    # Output sparse matrix with taxa_group/sample1/sample2/edge weight (number of supporting clusters with coverage > threshold)
    sparse_edges = sample_pairs.groupby([
        "taxa_group", "sample", "sample_2"
    ]).agg([
        pl.count().alias("weight").cast(int),
        pl.col("target").sort().str.concat(",").alias("target_ids")
    ]).select([
        "taxa_group", "weight", "target_ids",
        pl.col("sample").alias("sample1"),
        pl.col("sample_2").alias("sample2")
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

    unbinned = pl.read_csv(unbinned_path, sep="\t")

    targets, edges = pipeline(
        unbinned,
        MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
        TAXA_OF_INTEREST=TAXA_OF_INTEREST
        )
    targets.write_csv(targets_path, sep="\t")
    edges.write_csv(edges_path, sep="\t")
