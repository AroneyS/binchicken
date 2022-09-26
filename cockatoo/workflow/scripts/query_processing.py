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
for read in snakemake.input.reads:
    queries.append(pd.read_csv(read, sep = "\t"))

sequences = pd.concat(queries)

appraised = (sequences
    .reset_index()
    # Rename to match appraise output
    # Query output: query_name, query_sequence, divergence, num_hits, coverage, sample, marker, hit_sequence, taxonomy
    # Appraise output: gene, sample, sequence, num_hits, coverage, taxonomy, found_in
    .rename(columns = {"marker": "gene", "query_name":"sample", "query_sequence": "sequence", "sample": "found_in"})
    # Join query sample column to create found_in column
    .groupby(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "divergence"])["found_in"]
    .agg(lambda x: ",".join(sorted(x)))
    .reset_index()
    )

# Split dataframe into binned/unbinned
appraised["binned"] = appraised["divergence"].apply(lambda x: x <= (1 - SEQUENCE_IDENTITY) * WINDOW_SIZE)
binned = appraised[appraised["binned"]].drop(["divergence", "binned"], axis = 1)
unbinned = appraised[~appraised["binned"]].drop(["divergence", "binned"], axis = 1)

binned.to_csv(snakemake.output.binned, sep="\t", index=False)
unbinned.to_csv(snakemake.output.unbinned, sep="\t", index=False)
