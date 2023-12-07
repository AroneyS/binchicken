#############################
### precluster_samples.py ###
#############################
# Author: Samuel Aroney

import polars as pl
import os
import logging
from sourmash import MinHash
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

def processing(unbinned, MAX_CLUSTER_SIZE=1000):
    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    logging.info("Generating sketches")
    parent_mh = MinHash(n=0, ksize=60, scaled=1, track_abundance=False)
    hashes = []
    samples = []
    for group in unbinned.select("sample", "sequence").group_by("sample"):
        mh = parent_mh.copy_and_clear()
        for row in group[1].iter_rows():
            mh.add_sequence(row[1].replace("-", "A").replace("N", "A"))

        samples.append(group[0])
        hashes.append(mh)

    logging.info("Calculating distances")
    n_samples = len(hashes)
    distances = np.zeros([n_samples, n_samples])

    total_count = int(n_samples ** 2 / 2)
    for i, hash1 in enumerate(hashes):
        for j, hash2 in enumerate(hashes):
            if i < j:
                continue

            similarity = 1 - hash1.similarity(hash2)
            distances[i][j] = similarity
            distances[j][i] = similarity
    logging.info(f"Completed {total_count} comparisons")

    clust = linkage(squareform(distances), method="single")

    cluster_too_large = True
    t = 1
    while cluster_too_large:
        clusters = fcluster(clust, t=t)
        cluster_sizes = np.unique(clusters, return_counts=True)[1]
        cluster_too_large = np.any(cluster_sizes > MAX_CLUSTER_SIZE)
        t -= 0.1

    sample_clusters = []
    for cluster in np.unique(clusters):
        sample_cluster = [samples[s] for s in np.where(clusters == cluster)[0]]
        sample_clusters.append(sample_cluster)

    largest_cluster = max(sample_clusters, key=len)
    smallest_cluster = min(sample_clusters, key=len)
    logging.info(f"Found {len(sample_clusters)} clusters")
    logging.info(f"Largest cluster has {len(largest_cluster)} samples")
    logging.info(f"Smallest cluster has {len(smallest_cluster)} samples")

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

    unbinned_path = snakemake.input.unbinned
    KMER_PRECLUSTER = snakemake.params.kmer_precluster
    MAX_CLUSTER_SIZE = snakemake.params.max_precluster_size
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    output_dir = snakemake.output[0]
    os.makedirs(output_dir)

    unbinned = pl.read_csv(unbinned_path, separator="\t", dtypes=SINGLEM_OTU_TABLE_SCHEMA)

    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    if KMER_PRECLUSTER:
        clusters = processing(unbinned, MAX_CLUSTER_SIZE=MAX_CLUSTER_SIZE)

        with open(os.path.join(output_dir, "clusters.txt"), "w") as f:
            f.write("\n".join(sorted([",".join(sorted(cluster)) for cluster in clusters])) + "\n")

        for n, cluster in enumerate(clusters):
            (
                unbinned
                .filter(pl.col("sample").is_in(cluster))
                .write_csv(os.path.join(output_dir, "unbinned_" + str(n+1) + ".otu_table.tsv"), separator="\t")
            )
    else:
        unbinned.write_csv(os.path.join(output_dir, "unbinned_1.otu_table.tsv"), separator="\t")
