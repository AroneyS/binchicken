#!/usr/bin/env python3
###########################
### target_weighting.py ###
###########################
# Author: Samuel Aroney

import os
import polars as pl
import logging
import argparse

TARGET_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
    "target": str,
}

WEIGHTING_COLUMNS = {
    "gene": str,
    "sequence": str,
    "weight": float,
}

TARGET_WEIGHTING_COLUMNS = {
    "target": str,
    "weight": float,
}

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/target_weighting.py \
#   --targets {input.targets} \
#   --weighting {input.weighting} \
#   --targets-weighted {output.targets_weighted} \
#   --threads {threads} \
#   --log {log}
# """

def pipeline(targets, weighting):
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(targets) == 0:
        logging.warning("No targets found")
        return pl.DataFrame(schema=TARGET_WEIGHTING_COLUMNS)

    weighted = (
        targets
        .select("gene", "sequence", "target")
        .unique()
        .join(weighting, on=["gene", "sequence"])
        .drop("gene", "sequence")
    )

    return weighted

def main():
    parser = argparse.ArgumentParser(description="Target weighting pipeline script.")
    parser.add_argument("--targets", required=True, help="Path to input targets file")
    parser.add_argument("--weighting", required=True, help="Path to input weighting file")
    parser.add_argument("--targets-weighted", required=True, help="Path to output targets weighted file")
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

    targets = pl.read_csv(args.targets, separator="\t", schema_overrides=TARGET_COLUMNS)
    weighting = pl.read_csv(args.weighting, separator="\t", schema_overrides=WEIGHTING_COLUMNS)
    targets_weighted = pipeline(targets, weighting)
    targets_weighted.write_csv(args.targets_weighted, separator="\t")

if __name__ == "__main__":
    main()
