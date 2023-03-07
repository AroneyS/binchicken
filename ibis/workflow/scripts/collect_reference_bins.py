#################################
### collect_reference_bins.py ###
#################################
# Author: Samuel Aroney

import polars as pl
import os
import numpy as np
import extern

def trimmed_mean(data, trim=0.1):
    cut = int(np.floor(len(data) * trim))
    if cut == 0:
        return np.mean(data)
    else:
        a = sorted(data)
        return np.mean(a[cut:-cut])

def pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=0.1, TRIM_FRACTION=0.1):
    print(f"Polars using {str(pl.threadpool_size())} threads")

    appraise_binned = appraise_binned.with_columns(
        pl.col("sample").str.replace("\.1$", "")
    ).filter(
        pl.col("sample") == sample
    )

    appraise_unbinned = appraise_unbinned.with_columns(
        pl.col("sample").str.replace("\.1$", "")
    ).filter(
        pl.col("sample") == sample
    )

    num_binned = sum(appraise_binned.get_column("num_hits").to_list())
    num_unbinned = sum(appraise_unbinned.get_column("num_hits").to_list())
    perc_binned = num_binned / (num_binned + num_unbinned)

    if perc_binned < MIN_APPRAISED:
        return set()

    reference_bins = appraise_binned.with_columns(
        pl.col("found_in").str.split(",")
    ).explode(
        "found_in"
    ).with_columns(
        pl.col("found_in").str.replace("_protein$", "")
    ).groupby(
        ["gene", "found_in"]
    ).agg(
        pl.col("coverage").sum()
    ).pivot(
        values="coverage", index="gene", columns="found_in", 
    ).melt(
        id_vars="gene", variable_name="found_in", value_name="coverage"
    ).fill_null(0
    ).groupby(
        "found_in"
    ).agg(
        (pl.col("coverage").len() * TRIM_FRACTION).floor().cast(int).alias("cut"),
        pl.col("coverage")
    ).with_columns(
        pl.col("coverage").arr.sort().arr.slice(
            pl.col("cut"), pl.col("coverage").arr.lengths() - 2 * pl.col("cut")
            ).arr.mean()
    ).filter(
        pl.col("coverage") > 0
    ).get_column("found_in"
    ).to_list()

    return set(reference_bins)

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    binned_path = snakemake.input.appraise_binned
    unbinned_path = snakemake.input.appraise_unbinned
    genomes = snakemake.params.genomes
    sample = snakemake.params.sample
    MIN_APPRAISED = snakemake.params.min_appraised
    sample_read = snakemake.wildcards.read
    output_path = snakemake.output

    appraise_binned = pl.read_csv(binned_path, sep="\t")
    appraise_unbinned = pl.read_csv(unbinned_path, sep="\t")

    reference_bins = pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=MIN_APPRAISED)

    if len(reference_bins) == 0:
        print(f"Warning: No reference bins found for {sample_read}")
        cmd = f"touch {output_path}"
        extern.run(cmd)
    else:
        for bin in reference_bins:
            cmd = f"cat {genomes[bin]} >> {output_path}"
            extern.run(cmd)