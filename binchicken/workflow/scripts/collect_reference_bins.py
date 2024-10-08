#################################
### collect_reference_bins.py ###
#################################
# Author: Samuel Aroney

import polars as pl
import os
import numpy as np
import extern
from binchicken.binchicken import SUFFIX_RE

def trimmed_mean(data, trim=0.1):
    cut = int(np.floor(len(data) * trim))
    if cut == 0:
        return np.mean(data)
    else:
        a = sorted(data)
        return np.mean(a[cut:-cut])

def pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=0.1, TRIM_FRACTION=0.1):
    print(f"Polars using {str(pl.thread_pool_size())} threads")

    appraise_binned = (
        appraise_binned
        .with_columns(
            pl.col("sample").cast(str),
            sample_remove_suffix = pl.col("sample").str.replace(SUFFIX_RE, ""),
            )
        .filter((pl.col("sample") == sample) | (pl.col("sample_remove_suffix") == sample))
    )

    appraise_unbinned = (
        appraise_unbinned
        .with_columns(
            pl.col("sample").cast(str),
            sample_remove_suffix = pl.col("sample").str.replace(SUFFIX_RE, ""),
            )
        .filter((pl.col("sample") == sample) | (pl.col("sample_remove_suffix") == sample))
    )

    num_binned = sum(appraise_binned.get_column("num_hits").to_list())
    num_unbinned = sum(appraise_unbinned.get_column("num_hits").to_list())
    try:
        perc_binned = num_binned / (num_binned + num_unbinned)
    except ZeroDivisionError:
        return set()

    if perc_binned < MIN_APPRAISED:
        return set()

    reference_bins = (
        appraise_binned
        .with_columns(pl.col("found_in").str.split(","))
        .explode("found_in")
        .with_columns(pl.col("found_in").str.replace("_protein$", ""))
        .group_by(["gene", "found_in"])
        .agg(pl.col("coverage").sum())
        .pivot(values="coverage", index="gene", on="found_in", aggregate_function=None)
        .unpivot(index="gene", variable_name="found_in", value_name="coverage")
        .fill_null(0)
        .group_by("found_in")
        .agg(
            (pl.col("coverage").len() * TRIM_FRACTION).floor().cast(int).alias("cut"),
            pl.col("coverage")
            )
        .with_columns(
            pl.col("coverage")
                .list.sort()
                .list.slice(pl.col("cut"), pl.col("coverage").list.len() - 2 * pl.col("cut"))
                .list.mean()
            )
        .filter(pl.col("coverage") > 0)
        .get_column("found_in")
        .to_list()
    )

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

    appraise_binned = pl.read_csv(binned_path, separator="\t")
    appraise_unbinned = pl.read_csv(unbinned_path, separator="\t")

    reference_bins = pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=MIN_APPRAISED)

    if len(reference_bins) == 0:
        print(f"Warning: No reference bins found for {sample_read}")
        cmd = f"touch {output_path}"
        extern.run(cmd)
    else:
        with open(str(output_path), "w") as outfile:
            for bin in reference_bins:
                bin_name = os.path.splitext(os.path.basename(genomes[bin]))[0]
                with open(genomes[bin], "r") as infile:
                    for line in infile:
                        if line.startswith(">"):
                            line = f">{bin_name}~{line[1:]}"
                        outfile.write(line)
