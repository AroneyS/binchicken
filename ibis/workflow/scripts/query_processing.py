###########################
### query_processing.py ###
###########################
# Author: Samuel Aroney

import polars as pl
import os

OUTPUT_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": int,
    "taxonomy": str,
    "found_in": str,
    }

def processing(
    query_read,
    pipe_read,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60,
    TAXA_OF_INTEREST=None):

    if len(query_read) == 0:
        empty_output = pl.DataFrame(schema=OUTPUT_COLUMNS)
        return empty_output, empty_output

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        pipe_read = pipe_read.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    appraised = query_read.rename(
        # Rename to match appraise output
        # Query output: query_name, query_sequence, divergence, num_hits, coverage, sample, marker, hit_sequence, taxonomy
        # Appraise output: gene, sample, sequence, num_hits, coverage, taxonomy, found_in
        {"marker": "gene", "query_name":"sample", "query_sequence": "sequence", "sample": "found_in"}
    ).drop([
        # Query taxonomy, num_hits and coverage are from database (e.g. the genomes)
        "taxonomy", "num_hits", "coverage"
    ]).join(
        pipe_read, on=["gene", "sample", "sequence"], how="inner"
    ).groupby(
        ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "divergence"]
    ).agg(
        pl.col("found_in").sort().str.concat(",")
    ).with_columns(
        pl.col("divergence").alias("binned") <= ((1 - SEQUENCE_IDENTITY) * WINDOW_SIZE)
    )

    # Split dataframe into binned/unbinned
    binned = appraised.filter(pl.col("binned")).drop(["divergence", "binned"])
    unbinned = pipe_read.join(
        appraised, on=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"], how="left"
    ).filter(~pl.col("binned").fill_null(False)
    ).drop(["divergence", "binned"]
    ).with_columns(
        pl.lit(None).cast(str).alias("found_in")
    )

    # Filter out EIF
    binned = binned.filter(pl.col("gene") != "S3.18.EIF_2_alpha")
    unbinned = unbinned.filter(pl.col("gene") != "S3.18.EIF_2_alpha")

    return binned, unbinned

def pipeline(
    query_reads,
    pipe_reads,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60,
    TAXA_OF_INTEREST=None):

    print(f"Polars using {str(pl.threadpool_size())} threads")

    for query, pipe in zip(query_reads, pipe_reads):
        binned, unbinned = processing(
            pl.read_csv(query, separator="\t"),
            pl.read_csv(pipe, separator="\t"),
            SEQUENCE_IDENTITY,
            WINDOW_SIZE,
            TAXA_OF_INTEREST)
        yield binned, unbinned

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    SEQUENCE_IDENTITY = snakemake.params.sequence_identity
    WINDOW_SIZE = snakemake.params.window_size
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    query_reads = snakemake.input.query_reads
    pipe_reads = snakemake.input.pipe_reads
    binned_path = snakemake.output.binned
    unbinned_path = snakemake.output.unbinned

    outputs = pipeline(
        query_reads,
        pipe_reads,
        SEQUENCE_IDENTITY=SEQUENCE_IDENTITY,
        WINDOW_SIZE=WINDOW_SIZE,
        TAXA_OF_INTEREST=TAXA_OF_INTEREST
        )

    with open(binned_path, "ab") as binned_file, open(unbinned_path, "ab") as unbinned_file:
        first = True
        for binned, unbinned in outputs:
            binned.write_csv(binned_file, separator="\t", has_header=first)
            unbinned.write_csv(unbinned_file, separator="\t", has_header=first)
            first = False
