#!/usr/bin/env python3

"""
Author: Samuel Aroney
Check completion of Bin Chicken runs
"""

import os
import sys
import argparse
import logging
import polars as pl

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--dirs', nargs='+', help='Bin Chicken runs', required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    # Print whole string
    pl.Config(fmt_str_lengths=100)

    completion = (
        pl.DataFrame({"dir": args.dirs})
        .with_columns(
            elusive_clusters = pl.concat_str(pl.col("dir"), pl.lit("/coassemble/target/elusive_clusters.tsv")),
            )
        .with_columns(
            elusive_clusters = pl.col("elusive_clusters").map_elements(
                lambda x: pl.read_csv(x, separator="\t").select("coassembly").to_struct("elusive_clusters"),
                return_dtype = pl.List(pl.Struct([pl.Field("coassembly", pl.Utf8)]))
                ),
            )
        .explode("elusive_clusters")
        .unnest("elusive_clusters")
        .select("dir", "coassembly")
        .with_columns(
            assembled = pl.concat_str(pl.col("dir"), pl.lit("/coassemble/coassemble/"), pl.col("coassembly"), pl.lit("/assemble/assembly/final_contigs.fasta")),
            recovered = pl.concat_str(pl.col("dir"), pl.lit("/coassemble/coassemble/"), pl.col("coassembly"), pl.lit("/recover/bins/bin_info.tsv")),
            )
        .with_columns(
            assembled = pl.col("assembled").map_elements(lambda x: os.path.exists(x), return_dtype=pl.Boolean),
            recovered = pl.col("recovered").map_elements(lambda x: os.path.exists(x), return_dtype=pl.Boolean),
            )
    )

    summary = (
        completion
        .group_by("dir")
        .agg(
            total = pl.len(),
            assembled = pl.sum("assembled"),
            recovered = pl.sum("recovered"),
            )
        .with_columns(
            incomplete = pl.col("total") - pl.col("recovered"),
            )
    )

    toprint = (
        pl.concat([
            summary.sort("dir"),
            summary
                .group_by(len(args.dirs))
                .agg(
                    total = pl.sum("total"),
                    assembled = pl.sum("assembled"),
                    recovered = pl.sum("recovered"),
                    incomplete = pl.sum("incomplete"),
                    )
                .select(
                    pl.col("literal").cast(str).alias("dir"),
                    "total", "assembled", "recovered", "incomplete",
                    )
        ])
    )

    incomplete = completion.filter(pl.col("recovered").not_())

    logging.info("Incomplete:")
    for row in incomplete.iter_rows():
        try:
            completeness = "recover" if row[2] else "assemble"
            log_dir = os.path.join(row[0], "coassemble", "logs", "aviary", row[1] + "_" + completeness)
            log_subdir = os.path.join(log_dir, sorted(os.listdir(log_dir))[-1])
            log_path = sorted([os.path.join(log_subdir, f) for f in os.listdir(log_subdir) if f.endswith(".log")])[0]
            logging.info(f"{completeness} stage: {log_path}")
        except Exception as e:
            logging.info(f"Could not find log path for {row[0]} {row[1]}: {e}")

    logging.info(f"{toprint}")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
