#################################
### summarise_coassemblies.py ###
#################################
# Author: Samuel Aroney

import os
import polars as pl

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

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl
    
    elusive_clusters_path = snakemake.input.elusive_clusters
    read_size_path = snakemake.input.read_size
    output_path = snakemake.output.summary

    elusive_clusters = pl.read_csv(elusive_clusters_path, separator="\t")
    if read_size_path:
        read_size = pl.read_csv(read_size_path, has_header=False, new_columns=["sample", "read_size"])
    else:
        read_size = None

    summary = processing(elusive_clusters, read_size)
    summary.write_csv(output_path, separator="\t")
