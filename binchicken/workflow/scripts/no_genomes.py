#!/usr/bin/env python3
#####################
### no_genomes.py ###
#####################
# Author: Samuel Aroney

import polars as pl
import os
import argparse

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/no_genomes.py \
#   --reads {input.reads} \
#   --binned {output.binned} \
#   --unbinned {output.unbinned} \
#   --threads {threads}
# """

def processing(reads):
    print(f"Polars using {str(pl.thread_pool_size())} threads")

    unbinned = (
        reads
        .with_columns(found_in = pl.lit(""))
        .filter(pl.col("gene") != "S3.18.EIF_2_alpha")
    )
    binned = unbinned.filter(False)

    return(binned, unbinned)

def main():
    parser = argparse.ArgumentParser(description="No genomes pipeline script.")
    parser.add_argument("--reads", required=True, help="List file of input read files")
    parser.add_argument("--binned", required=True, help="Path to output binned file")
    parser.add_argument("--unbinned", required=True, help="Path to output unbinned file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    reads = pl.scan_csv(args.reads, separator="\t", schema_overrides=SINGLEM_OTU_TABLE_SCHEMA)
    binned, unbinned = processing(reads)

    binned.sink_csv(args.binned, separator="\t")
    unbinned.sink_csv(args.unbinned, separator="\t")

if __name__ == "__main__":
    main()
