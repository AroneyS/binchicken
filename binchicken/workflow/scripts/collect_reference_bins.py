#!/usr/bin/env python3
#################################
### collect_reference_bins.py ###
#################################
# Author: Samuel Aroney

import polars as pl
import os
import numpy as np
import extern
import argparse

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/collect_reference_bins.py \
#   --appraise-binned {input.appraise_binned} \
#   --appraise-unbinned {input.appraise_unbinned} \
#   --genomes {input.genomes} \
#   --sample {params.sample} \
#   --min-appraised {params.min_appraised} \
#   --read {wildcards.read} \
#   --output {output}
# """

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
        .with_columns(pl.col("sample").cast(str))
        .filter(pl.col("sample") == sample)
    )

    appraise_unbinned = (
        appraise_unbinned
        .with_columns(pl.col("sample").cast(str))
        .filter(pl.col("sample") == sample)
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

def main():
    parser = argparse.ArgumentParser(description="Collect reference bins pipeline script.")
    parser.add_argument("--appraise-binned", required=True, help="Path to appraise binned input file")
    parser.add_argument("--appraise-unbinned", required=True, help="Path to appraise unbinned input file")
    parser.add_argument("--genomes", required=True, help="Named list file of genome fasta files (indexed by bin name)")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--min-appraised", type=float, default=0.1, help="Minimum appraised fraction")
    parser.add_argument("--read", required=True, help="Read wildcard value")
    parser.add_argument("--output", required=True, help="Output fasta file path")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    appraise_binned = pl.read_csv(args.appraise_binned, separator="\t")
    appraise_unbinned = pl.read_csv(args.appraise_unbinned, separator="\t")

    # Build a mapping from bin name to genome file path
    genomes_dict = {}
    with open(args.genomes, "r") as f:
        for line in f:
            bin_name, genome_path = line.strip().split("\t")
            genomes_dict[bin_name] = genome_path

    reference_bins = pipeline(appraise_binned, appraise_unbinned, args.sample, MIN_APPRAISED=args.min_appraised)

    if len(reference_bins) == 0:
        print(f"Warning: No reference bins found for {args.read}")
        cmd = f"touch {args.output}"
        extern.run(cmd)
    else:
        with open(str(args.output), "w") as outfile:
            for bin in reference_bins:
                bin_name = os.path.splitext(os.path.basename(genomes_dict[bin]))[0]
                with open(genomes_dict[bin], "r") as infile:
                    for line in infile:
                        if line.startswith(">"):
                            line = f">{bin_name}~{line[1:]}"
                        outfile.write(line)

if __name__ == "__main__":
    main()
