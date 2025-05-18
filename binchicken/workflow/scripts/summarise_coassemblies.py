#!/usr/bin/env python3
#################################
### summarise_coassemblies.py ###
#################################
# Author: Samuel Aroney

import os
import polars as pl
import argparse

def processing(elusive_clusters, read_size):
    print(f"Polars using {str(pl.thread_pool_size())} threads")

    summary = (
        elusive_clusters
        .select("coassembly", "samples", "length", "total_targets", "total_size")
    )

    if read_size is not None:
        summary = (
            summary
            .with_columns(sample = pl.col("samples").str.split(","))
            .explode("sample")
            .join(read_size, on="sample", how="left", coalesce=True)
            .group_by("coassembly", "samples", "length", "total_targets", "total_size")
            .agg(unmapped_size = pl.sum("read_size"))
        )

    return summary

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/summarise_coassemblies.py \
#   --elusive-clusters {input.elusive_clusters} \
#   --read-size {input.read_size} \
#   --summary {output.summary} \
#   --threads {threads}
# """

def main():
    parser = argparse.ArgumentParser(description="Summarise coassemblies pipeline script.")
    parser.add_argument("--elusive-clusters", required=True, help="Path to elusive clusters input file")
    parser.add_argument("--read-size", nargs='?', const=None, default=None, help="Path to read size input file (optional)")
    parser.add_argument("--summary", required=True, help="Path to output summary file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    elusive_clusters = pl.read_csv(args.elusive_clusters, separator="\t")
    if args.read_size:
        read_size = pl.read_csv(args.read_size, has_header=False, new_columns=["sample", "read_size"])
    else:
        read_size = None

    summary = processing(elusive_clusters, read_size)
    summary.write_csv(args.summary, separator="\t")

if __name__ == "__main__":
    main()
