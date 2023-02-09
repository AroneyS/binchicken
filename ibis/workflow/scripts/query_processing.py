###########################
### query_processing.py ###
###########################
# Author: Samuel Aroney

import pandas as pd

def processing(
    query_read,
    pipe_read,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60):

    output_columns = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

    if len(query_read) == 0:
        empty_output = pd.DataFrame(columns=output_columns)
        return empty_output, empty_output

    appraised = (query_read
        # Rename to match appraise output
        # Query output: query_name, query_sequence, divergence, num_hits, coverage, sample, marker, hit_sequence, taxonomy
        # Appraise output: gene, sample, sequence, num_hits, coverage, taxonomy, found_in
        .rename(columns = {"marker": "gene", "query_name":"sample", "query_sequence": "sequence", "sample": "found_in"})
        # Query taxonomy, num_hits and coverage are from database (e.g. the genomes)
        .drop(["taxonomy", "num_hits", "coverage"], axis=1, errors="ignore")
        .set_index(["gene", "sample", "sequence"])
        .join(pipe_read.set_index(["gene", "sample", "sequence"]), on = ["gene", "sample", "sequence"], how = "inner")
        .reset_index()
        # Join query sample column to create found_in column
        .groupby(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "divergence"])["found_in"]
        .agg(lambda x: ",".join(sorted(x)))
        .reset_index()
        )

    # Split dataframe into binned/unbinned
    appraised["binned"] = appraised["divergence"].apply(lambda x: x <= (1 - SEQUENCE_IDENTITY) * WINDOW_SIZE)
    binned = appraised[appraised["binned"]].drop(["divergence", "binned"], axis = 1).reset_index(drop = True)
    unbinned = pipe_read.join(appraised.set_index(output_columns[0:-1]), on = output_columns[0:-1])
    unbinned = unbinned[~unbinned["binned"].fillna(False)].drop(["divergence", "binned"], axis = 1).reset_index(drop = True)
    unbinned["found_in"] = None

    return binned, unbinned

def pipeline(
    query_reads,
    pipe_reads,
    SEQUENCE_IDENTITY=0.86,
    WINDOW_SIZE=60):

    for query, pipe in zip(query_reads, pipe_reads):
        binned, unbinned = processing(
            pd.read_csv(query, sep = "\t"),
            pd.read_csv(pipe, sep = "\t"),
            SEQUENCE_IDENTITY,
            WINDOW_SIZE)
        yield binned, unbinned

if __name__ == "__main__":
    SEQUENCE_IDENTITY = snakemake.params.sequence_identity
    WINDOW_SIZE = snakemake.params.window_size
    query_reads = snakemake.input.query_reads
    pipe_reads = snakemake.input.pipe_reads
    binned_path = snakemake.output.binned
    unbinned_path = snakemake.output.unbinned

    outputs = pipeline(
        query_reads,
        pipe_reads,
        SEQUENCE_IDENTITY=SEQUENCE_IDENTITY,
        WINDOW_SIZE=WINDOW_SIZE
        )

    first = True
    for binned, unbinned in outputs:
        binned.to_csv(binned_path, sep = "\t", mode = "a", header = first, index = False)
        unbinned.to_csv(unbinned_path, sep = "\t", mode = "a", header = first, index = False)
        first = False
