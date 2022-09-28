###########################
### query_processing.py ###
###########################
# Author: Samuel Aroney

import pandas as pd

SEQUENCE_IDENTITY = snakemake.params.sequence_identity
WINDOW_SIZE = snakemake.params.window_size

################
### Pipeline ###
################
queries = []
for read in snakemake.input.query_reads:
    queries.append(pd.read_csv(read, sep = "\t"))
query_sequences = pd.concat(queries)

pipes = []
for read in snakemake.input.pipe_reads:
    pipes.append(pd.read_csv(read, sep = "\t"))
pipe_sequences = pd.concat(pipes)

appraised = (query_sequences
    # Rename to match appraise output
    # Query output: query_name, query_sequence, divergence, num_hits, coverage, sample, marker, hit_sequence, taxonomy
    # Appraise output: gene, sample, sequence, num_hits, coverage, taxonomy, found_in
    .rename(columns = {"marker": "gene", "query_name":"sample", "query_sequence": "sequence", "sample": "found_in"})
    # Query taxonomy, num_hits and coverage are from database (e.g. the genomes)
    .drop(["taxonomy", "num_hits", "coverage"], axis = 1)
    .set_index(["gene", "sample", "sequence"])
    .join(pipe_sequences.set_index(["gene", "sample", "sequence"]), on = ["gene", "sample", "sequence"])
    .reset_index()
    # Join query sample column to create found_in column
    .groupby(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "divergence"])["found_in"]
    .agg(lambda x: ",".join(sorted(x)))
    .reset_index()
    )

# Split dataframe into binned/unbinned
appraised["binned"] = appraised["divergence"].apply(lambda x: x <= (1 - SEQUENCE_IDENTITY) * WINDOW_SIZE)
binned = appraised[appraised["binned"]].drop(["divergence", "binned"], axis = 1)
unbinned = appraised[~appraised["binned"]].drop(["divergence", "binned"], axis = 1)
unbinned["found_in"] = None

binned.to_csv(snakemake.output.binned, sep="\t", index=False)
unbinned.to_csv(snakemake.output.unbinned, sep="\t", index=False)
