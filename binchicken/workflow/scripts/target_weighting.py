###########################
### target_weighting.py ###
###########################
# Author: Samuel Aroney

import os
import polars as pl
import logging

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

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    targets_path = snakemake.input.targets
    weighting_path = snakemake.input.weighting
    targets_weighted_path = snakemake.output.targets_weighted

    targets = pl.read_csv(targets_path, separator="\t", dtypes=TARGET_COLUMNS)
    weighting = pl.read_csv(weighting_path, separator="\t", dtypes=WEIGHTING_COLUMNS)

    targets_weighted = pipeline(targets, weighting)
    targets_weighted.write_csv(targets_weighted_path, separator="\t")
