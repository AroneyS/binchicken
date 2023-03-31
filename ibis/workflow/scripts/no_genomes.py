#####################
### no_genomes.py ###
#####################
# Author: Samuel Aroney

import pandas as pd

reads = []
for read in snakemake.input.reads:
    reads.append(pd.read_csv(read, sep = "\t"))

unbinned = pd.concat(reads)
unbinned["found_in"] = ""
unbinned.to_csv(snakemake.output.unbinned, sep="\t", index=False)

binned = unbinned.drop(unbinned.index)
binned.to_csv(snakemake.output.binned, sep="\t", index=False)
