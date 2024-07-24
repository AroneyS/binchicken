##############################
### abundance_weighting.py ###
##############################
# Author: Samuel Aroney

import os
import polars as pl
import logging

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

    weighted = (
        unbinned
        .join(total_coverage, on=["sample", "gene"])
        .with_columns(weight = pl.col("coverage") / pl.col("total_coverage"))
        .group_by("gene", "sequence")
        .agg(pl.mean("weight"))
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

    unbinned_path = snakemake.input.unbinned
    binned_path = snakemake.input.binned
    weighted_path = snakemake.output.weighted
    samples = snakemake.params.samples

    unbinned = pl.read_csv(unbinned_path, separator="\t", dtypes=APPRAISE_COLUMNS)
    binned = pl.read_csv(binned_path, separator="\t", dtypes=APPRAISE_COLUMNS)

    weighted = pipeline(unbinned, binned, samples)
    weighted.write_csv(weighted_path, separator="\t")
