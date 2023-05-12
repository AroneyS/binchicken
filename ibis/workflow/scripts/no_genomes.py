#####################
### no_genomes.py ###
#####################
# Author: Samuel Aroney

import polars as pl
import os

def processing(reads):
    print(f"Polars using {str(pl.threadpool_size())} threads")

    unbinned = (
        reads
        .with_columns(found_in = pl.lit(""))
        .filter(pl.col("gene") != "S3.18.EIF_2_alpha")
    )
    binned = unbinned.filter(False)

    return(binned, unbinned)

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    input_reads = snakemake.input.reads
    output_binned = snakemake.output.binned
    output_unbinned = snakemake.output.unbinned

    reads = []
    for read in input_reads:
        reads.append(pl.read_csv(read, separator="\t"))

    binned, unbinned = processing(pl.concat(reads))

    binned.write_csv(output_binned, separator="\t")
    unbinned.write_csv(output_unbinned, separator="\t")
