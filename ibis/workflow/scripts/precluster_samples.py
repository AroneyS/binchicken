#############################
### precluster_samples.py ###
#############################
# Author: Samuel Aroney

import polars as pl
import os
import subprocess

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

def processing(unbinned, working_dir):
    print(f"Polars using {str(pl.threadpool_size())} threads")
    os.makedirs(working_dir)

    # Generate fa files
    for group in unbinned.select("sample", "sequence").group_by("sample"):
        sample_name = group[0]
        with open(os.path.join(working_dir, sample_name + ".fa"), "w") as f:
            for n, row in enumerate(group[1].iter_rows()):
                if n == 0:
                    seq_name = sample_name
                else:
                    seq_name = sample_name + str(n)

                seq = row[1].replace("-", "A")
                f.write(f">{seq_name}\n{seq}\n")

    with open(os.path.join(working_dir, "samples.txt"), "w") as f:
        for sample in unbinned.unique("sample").get_column("sample"):
            filename = os.path.join(working_dir, sample + ".fa")
            f.write(f"{filename}\n")

    # Generate large kmers
    cmd = (
        f"sourmash sketch dna "
        f"-p k=60,scaled=1,noabund "
        f"--name-from-first "
        f"--from-file {working_dir}/samples.txt "
        f"-o {working_dir}/sketches.sig "
    )
    subprocess.run(cmd, shell=True, check=True)

    cmd = (
        f"sourmash compare {working_dir}/sketches.sig -o {working_dir}/sourmash_compare && "
        f"sourmash plot --labels {working_dir}/sourmash_compare"
    )
    subprocess.run(cmd, shell=True, check=True)

    # Cluster into groups of max 1000

    return(clusters)

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    unbinned_path = snakemake.input.unbinned
    output_dir = snakemake.output[0]

    unbinned = pl.scan_csv(unbinned_path, separator="\t", dtypes=SINGLEM_OTU_TABLE_SCHEMA)
    clusters = processing(unbinned, os.path.join(output_dir, "working_dir"))

    for n, cluster in enumerate(clusters.iter_rows()):
        (
            unbinned
            .filter(pl.col("sample").isin(cluster))
            .write_csv(os.path.join(output_dir, "unbinned_" + str(n+1) + ".otu_table.tsv"), separator="\t")
        )
