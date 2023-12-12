#############################
### precluster_samples.py ###
#############################
# Author: Samuel Aroney

import polars as pl
import os
import shutil
import logging
from sourmash import fig
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

def processing(distances, samples, MAX_CLUSTER_SIZE=1000):
    logging.info(f"Clustering samples")
    clust = linkage(squareform(distances), method="single",  metric="precomputed")

    cluster_too_large = True
    t_increment = 0.01
    t = 1 + t_increment
    while cluster_too_large:
        t -= t_increment
        clusters = fcluster(clust, criterion="distance", t=t)
        cluster_sizes = np.unique(clusters, return_counts=True)[1]
        cluster_too_large = np.any(cluster_sizes > MAX_CLUSTER_SIZE)

    logging.info(f"Found cutoff t={round(t, ndigits=2)} with no clusters larger than {MAX_CLUSTER_SIZE}")

    sample_clusters = []
    for cluster in np.unique(clusters):
        sample_cluster = [samples[s] for s in np.where(clusters == cluster)[0]]
        sample_clusters.append(sample_cluster)

    logging.info(f"Found {len(sample_clusters)} clusters")
    logging.info(f"Largest cluster has {max(cluster_sizes)} samples")
    logging.info(f"Smallest cluster has {min(cluster_sizes)} samples")

    return sample_clusters

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    distances_path = snakemake.input.distance
    unbinned_path = snakemake.input.unbinned
    KMER_PRECLUSTER = snakemake.params.kmer_precluster
    MAX_CLUSTER_SIZE = snakemake.params.max_precluster_size
    output_dir = snakemake.output[0]
    os.makedirs(output_dir)

    if KMER_PRECLUSTER:
        distances, samples = fig.load_matrix_and_labels(distances_path)
        clusters = processing(distances, samples, MAX_CLUSTER_SIZE=MAX_CLUSTER_SIZE)

        with open(os.path.join(output_dir, "clusters.txt"), "w") as f:
            f.write("\n".join(sorted([",".join(sorted(cluster)) for cluster in clusters])) + "\n")

        unbinned = pl.read_csv(unbinned_path, separator="\t", dtypes=SINGLEM_OTU_TABLE_SCHEMA)
        for n, cluster in enumerate(clusters):
            (
                unbinned
                .filter(pl.col("sample").is_in(cluster))
                .write_csv(os.path.join(output_dir, "unbinned_" + str(n+1) + ".otu_table.tsv"), separator="\t")
            )
    else:
        shutil.copyfile(unbinned_path, os.path.join(output_dir, "unbinned_1.otu_table.tsv"))
