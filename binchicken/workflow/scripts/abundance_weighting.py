#!/usr/bin/env python3
##############################
### abundance_weighting.py ###
##############################
# Author: Samuel Aroney

import os
import polars as pl
import logging
import argparse

APPRAISE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
}

WEIGHTING_COLUMNS = {
    "gene": str,
    "sequence": str,
    "weight": float,
}

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/abundance_weighting.py \
#   --unbinned {input.unbinned} \
#   --binned {input.binned} \
#   --weighted {output.weighted} \
#   --samples {input.samples} \
#   --threads {threads} \
#   --log {log}
# """

def pipeline(unbinned, binned, samples=None):
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(unbinned) == 0:
        logging.warning("No unbinned sequences found")
        return pl.DataFrame(schema=WEIGHTING_COLUMNS)

    total_coverage = (
        pl.concat([binned, unbinned])
        .group_by("sample", "gene")
        .agg(total_coverage = pl.sum("coverage"))
    )

    if samples:
        total_coverage = (
            total_coverage
            .filter(pl.col("sample").is_in(samples))
        )

    gene_sequences = (
        unbinned
        .select("gene", "sequence")
        .unique()
    )

    weighted = (
        total_coverage
        .join(gene_sequences, on=["gene"])
        .join(unbinned, on=["sample", "gene", "sequence"], how="full")
        .filter(pl.col("sample").is_not_null())
        .with_columns(pl.col("coverage").fill_null(0))
        .with_columns(weight = pl.col("coverage") / pl.col("total_coverage"))
        .group_by("gene", "sequence")
        .agg(pl.mean("weight"))
        .filter(pl.col("weight") > 0)
    )

    return weighted

def main():
    parser = argparse.ArgumentParser(description="Abundance weighting pipeline script.")
    parser.add_argument("--unbinned", required=True, help="Path to unbinned input file")
    parser.add_argument("--binned", required=True, help="Path to binned input file")
    parser.add_argument("--weighted", required=True, help="Path to output weighted file")
    parser.add_argument("--samples", default=None, help="List file of sample names")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if args.log:
        logging.basicConfig(
            filename=args.log,
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    unbinned = pl.read_csv(args.unbinned, separator="\t", schema_overrides=APPRAISE_COLUMNS)
    binned = pl.read_csv(args.binned, separator="\t", schema_overrides=APPRAISE_COLUMNS)

    if not args.samples:
        samples = None
    elif os.path.getsize(args.samples) == 0:
        samples = None
    else:
        samples = pl.read_csv(args.samples, has_header=False, new_columns=["sample"]).get_column("sample").to_list()

    weighted = pipeline(unbinned, binned, samples)
    weighted.write_csv(args.weighted, separator="\t")

if __name__ == "__main__":
    main()
