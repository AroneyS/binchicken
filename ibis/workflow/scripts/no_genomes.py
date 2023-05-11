#####################
### no_genomes.py ###
#####################
# Author: Samuel Aroney

import polars as pl

reads = []
for read in snakemake.input.reads:
    reads.append(pl.read_csv(read, separator="\t"))

unbinned = pl.concat(reads).with_columns(found_in = pl.lit(""))
unbinned.write_csv(snakemake.output.unbinned, separator="\t")

binned = unbinned.filter(False)
binned.write_csv(snakemake.output.binned, separator="\t")
