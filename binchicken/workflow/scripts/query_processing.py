#!/usr/bin/env python3
###########################
### query_processing.py ###
###########################
# Author: Samuel Aroney

import polars as pl
import os
import logging
import argparse

QUERY_COLUMNS = {
    "query_name": str,
    "query_sequence": str,
    "divergence": int,
    "num_hits": int,
    "coverage": float,
    "sample": str,
    "marker": str,
    "hit_sequence": str,
    "taxonomy": str,
}

PIPE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": int,
    "taxonomy": str,
}

OUTPUT_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": int,
    "taxonomy": str,
    "found_in": str,
}

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/query_processing.py \
#   --query-reads {input.query_reads} \
#   --pipe-reads {input.pipe_reads} \
#   --binned {output.binned} \
#   --unbinned {output.unbinned} \
#   --sequence-identity {params.sequence_identity} \
#   --window-size {params.window_size} \
#   --threads {threads} \
#   --log {log}
# """

def processing(
    query_read,
    pipe_read,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60):

    if len(query_read) == 0:
        empty_output = pl.DataFrame(schema=OUTPUT_COLUMNS)
        return empty_output, empty_output

    appraised = (
        query_read
        # Rename to match appraise output
        # Query output: query_name, query_sequence, divergence, num_hits, coverage, sample, marker, hit_sequence, taxonomy
        # Appraise output: gene, sample, sequence, num_hits, coverage, taxonomy, found_in
        .rename({"marker": "gene", "query_name":"sample", "query_sequence": "sequence", "sample": "found_in"})
        # Query taxonomy, num_hits and coverage are from database (i.e. the reference genomes)
        .drop(["taxonomy", "num_hits", "coverage"])
        .join(pipe_read, on=["gene", "sample", "sequence"], how="inner")
        .group_by(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "divergence"])
        .agg(pl.col("found_in").sort().str.concat(","))
        .with_columns(pl.col("divergence").alias("binned") <= ((1 - SEQUENCE_IDENTITY) * WINDOW_SIZE))
    )

    # Split dataframe into binned/unbinned
    binned = appraised.filter(pl.col("binned")).drop(["divergence", "binned"])
    unbinned = (
        pipe_read
        .join(
            appraised, on=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"], how="left", coalesce=True
            )
        .filter(~pl.col("binned").fill_null(False))
        .drop(["divergence", "binned"])
        .with_columns(pl.lit(None).cast(str).alias("found_in"))
    )

    return binned, unbinned

def pipeline(
    query_reads,
    pipe_reads,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60):

    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    for query, pipe in zip(query_reads, pipe_reads):
        logging.debug(f"Processing {query} and {pipe}")
        binned, unbinned = processing(
            pl.read_csv(query, separator="\t"),
            pl.read_csv(pipe, separator="\t"),
            SEQUENCE_IDENTITY,
            WINDOW_SIZE)
        yield binned, unbinned

def main():
    parser = argparse.ArgumentParser(description="Query processing pipeline script.")
    parser.add_argument("--query-reads", required=True, nargs='+', help="List of query read files")
    parser.add_argument("--pipe-reads", required=True, nargs='+', help="List of pipe read files")
    parser.add_argument("--binned", required=True, help="Path to output binned file")
    parser.add_argument("--unbinned", required=True, help="Path to output unbinned file")
    parser.add_argument("--sequence-identity", type=float, default=0.86, help="Sequence identity threshold")
    parser.add_argument("--window-size", type=int, default=60, help="Window size")
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

    outputs = pipeline(
        args.query_reads,
        args.pipe_reads,
        SEQUENCE_IDENTITY=args.sequence_identity,
        WINDOW_SIZE=args.window_size
    )

    with open(args.binned, "ab") as binned_file, open(args.unbinned, "ab") as unbinned_file:
        first = True
        for binned, unbinned in outputs:
            binned.write_csv(binned_file, separator="\t", include_header=first)
            unbinned.write_csv(unbinned_file, separator="\t", include_header=first)
            first = False

    logging.info("Done")

if __name__ == "__main__":
    main()
