#########################
### target_elusive.py ###
#########################
# Author: Samuel Aroney

import os
import polars as pl
import logging
import numpy as np
import scipy.sparse as sp
import itertools

EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
    }

def get_clusters(
        sample_distances,
        PRECLUSTER_SIZE=2,
        MAX_COASSEMBLY_SAMPLES=2):
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if sample_distances.height == 0:
        return pl.DataFrame(schema={"samples": str})

    if MAX_COASSEMBLY_SAMPLES < 2:
        # Set to 2 to produce paired edges
        MAX_COASSEMBLY_SAMPLES = 2

    logging.info("Converting to sparse array")
    samples = np.unique(np.concatenate([
        sample_distances.select("query_name").to_numpy().flatten(),
        sample_distances.select("match_name").to_numpy().flatten(),
    ]))
    sample_to_index = {sample: index for index, sample in enumerate(samples)}

    logging.info("Initialise the array")
    distances = (
        sp.coo_matrix(
            (
                sample_distances.select("jaccard").to_numpy().flatten().astype(np.float32),
                (
                    np.array([sample_to_index[name] for name in sample_distances.select("query_name").to_numpy().flatten()]),
                    np.array([sample_to_index[name] for name in sample_distances.select("match_name").to_numpy().flatten()])
                )
            ),
            shape=(len(samples), len(samples))
        )
        .tocsr()
    )
    distances = distances + distances.transpose()
    distances.setdiag(1)

    logging.info("Processing distances...")
    best_samples = np.empty((distances.shape[0], PRECLUSTER_SIZE), dtype=np.int32)
    for i in range(distances.shape[0]):
        row = distances.getrow(i).toarray().ravel()
        best_samples[i] = np.argsort(row)[-PRECLUSTER_SIZE:]

    chosen_samples = [[[samples[i]] + list(np.array(samples)[b[b != i]])] for i, b in enumerate(best_samples)]

    sample_combinations = (
        pl.LazyFrame({"cluster_size": range(1, MAX_COASSEMBLY_SAMPLES)})
        # Choose top N clusters for each cluster size where N is PRECLUSTER_SIZE
        .with_columns(
            sample_combinations = pl.col("cluster_size").map_elements(
                lambda x: [[0] + list(i) for i in itertools.combinations(range(1, PRECLUSTER_SIZE), x)][:PRECLUSTER_SIZE],
                return_dtype=pl.List(pl.List(pl.Int64)),
                )
            )
        .explode("sample_combinations")
        .select("sample_combinations")
    )

    logging.info("Choosing preclusters based on distances")
    with pl.StringCache():
        preclusters = (
            pl.LazyFrame(chosen_samples, schema={"samples": pl.List(pl.Categorical)}, orient="row")
            .join(sample_combinations, how="cross")
            .select(
                pl.col("samples")
                    .list.gather(pl.col("sample_combinations"))
                    .cast(pl.List(str))
                    .list.sort()
                    .list.join(",")
                )
            .unique()
            .collect(streaming=True)
        )

    logging.info(f"Found {preclusters.height} preclusters")
    return preclusters

def streaming_pipeline(
    unbinned,
    samples,
    sample_preclusters,
    targets_path,
    edges_path,
    MIN_COASSEMBLY_COVERAGE=10,
    TAXA_OF_INTEREST="",
    MAX_COASSEMBLY_SAMPLES=2,
    CHUNK_SIZE=2):

    # logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(unbinned) == 0:
        logging.warning("No unbinned sequences found")
        unbinned.rename({"found_in": "target"}).write_csv(targets_path, separator="\t")
        pl.DataFrame(schema=EDGES_COLUMNS).write_csv(edges_path, separator="\t")
        return

    if MAX_COASSEMBLY_SAMPLES < 2:
        # Set to 2 to produce paired edges
        MAX_COASSEMBLY_SAMPLES = 2

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    logging.info("Grouping hits by marker gene sequences to form targets")
    unbinned = (
        unbinned
        .with_columns(
            pl.when(pl.col("sample").is_in(samples))
            .then(pl.col("sample"))
            .otherwise(pl.col("sample").str.replace(r"(_|\.)R?1$", ""))
            )
        .filter(pl.col("sample").is_in(samples))
        .drop("found_in")
        .with_row_index("target")
        .select(
            "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
            pl.first("target").over(["gene", "sequence"]).rank("dense") - 1,
            )
    )

    unbinned.with_columns(pl.col("target").cast(pl.Utf8)).write_csv(targets_path, separator="\t")

    if unbinned.height == 0:
        logging.warning("No SingleM sequences found for the given samples")
        pl.DataFrame(schema=EDGES_COLUMNS).write_csv(edges_path, separator="\t")
        return

    if sample_preclusters.height == 0:
        logging.warning("No preclusters found")
        pl.DataFrame(schema=EDGES_COLUMNS).write_csv(edges_path, separator="\t")
        return

    logging.info("Using chosen clusters to find appropriate targets")
    def process_chunk(df):
        sparse_edges = (
            df
            .lazy()
            .with_columns(sample_ids = pl.col("samples").str.split(",").cast(pl.List(pl.Categorical)))
            .with_columns(cluster_size = pl.col("sample_ids").list.len())
            .explode("sample_ids")
            .join(unbinned.lazy().select("target", "coverage", sample_ids=pl.col("sample").cast(pl.Categorical)), on="sample_ids")
            .group_by("samples", "cluster_size", "target")
            .agg(pl.sum("coverage"), count = pl.len())
            .filter(pl.col("count") == pl.col("cluster_size"))
            .filter(pl.col("coverage") > MIN_COASSEMBLY_COVERAGE)
            .group_by("samples", "cluster_size")
            .agg(target_ids = pl.col("target").cast(pl.Utf8).sort().str.concat(","))
            .with_columns(style = pl.lit("match"))
            .select("style", "cluster_size", "samples", "target_ids")
            .collect(streaming=True)
        )

        return(sparse_edges)

    num_chunks = (sample_preclusters.height + CHUNK_SIZE - 1) // CHUNK_SIZE # Ceiling division to include all rows
    with open(edges_path, "w") as f:
        with pl.StringCache():
            for i in range(num_chunks):
                start_row = i * CHUNK_SIZE
                chunk = sample_preclusters.slice(start_row, CHUNK_SIZE)
                processed_chunk = process_chunk(chunk)
                processed_chunk.write_csv(f, separator="\t", include_header=i==0)

    return

def pipeline(
    unbinned,
    samples,
    MIN_COASSEMBLY_COVERAGE=10,
    TAXA_OF_INTEREST="",
    MAX_COASSEMBLY_SAMPLES=2):

    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    if len(unbinned) == 0:
        logging.warning("No unbinned sequences found")
        return unbinned.rename({"found_in": "target"}), pl.DataFrame(schema=EDGES_COLUMNS)

    if MAX_COASSEMBLY_SAMPLES < 2:
        # Set to 2 to produce paired edges
        MAX_COASSEMBLY_SAMPLES = 2

    # Filter TAXA_OF_INTEREST
    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    logging.info("Grouping hits by marker gene sequences to form targets")
    unbinned = (
        unbinned
        .with_columns(
            pl.when(pl.col("sample").is_in(samples))
            .then(pl.col("sample"))
            .otherwise(pl.col("sample").str.replace(r"(_|\.)R?1$", ""))
            )
        .filter(pl.col("sample").is_in(samples))
        .drop("found_in")
        .with_row_index("target")
        .select(
            "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
            pl.first("target").over(["gene", "sequence"]).rank("dense") - 1,
            )
    )

    if unbinned.height == 0:
        logging.warning("No SingleM sequences found for the given samples")
        return unbinned.with_columns(pl.col("target").cast(pl.Utf8)), pl.DataFrame(schema=EDGES_COLUMNS)

    def process_groups(df):
        if df.height == 1:
            return pl.DataFrame(schema={"style": str, "cluster_size": pl.Int64, "samples": str, "target": pl.UInt32})

        # Direct matching samples in pairs with coverage > MIN_COASSEMBLY_COVERAGE
        dfs = [
            df
            .join(df, how="cross", suffix="_2")
            .filter(
                pl.col("samples").list.last().str.encode("hex") < pl.col("samples_2").list.first().str.encode("hex")
                )
            .select([
                "target",
                pl.col("coverage") + pl.col("coverage_2").alias("coverage"),
                pl.col("samples").list.concat(pl.col("samples_2")),
                ])
            .filter(pl.col("samples").list.len() > 1)
            .filter(pl.col("coverage") > MIN_COASSEMBLY_COVERAGE)
            .select(
                pl.lit("match").alias("style"),
                pl.lit(2).cast(pl.Int64).alias("cluster_size"),
                pl.col("samples").list.join(","),
                pl.col("target"),
                )
        ]

        # Pool samples with coverage > MIN_COASSEMBLY_COVERAGE / N for clusters of size N: 3 to MAX_COASSEMBLY_SAMPLES
        cluster_sizes = pl.DataFrame({"cluster_size": range(3, MAX_COASSEMBLY_SAMPLES+1)})
        dfs.append(
            df
            .join(cluster_sizes, how="cross")
            .filter(pl.col("coverage") > float(MIN_COASSEMBLY_COVERAGE) / pl.col("cluster_size").cast(float))
            .group_by("target", "cluster_size")
            .agg(pl.concat_list("samples").flatten())
            .filter(pl.col("samples").list.len() >= pl.col("cluster_size"))
            .select(
                pl.lit("pool").alias("style"),
                pl.col("cluster_size"),
                pl.col("samples").list.sort().list.join(","),
                pl.col("target"),
                )
        )

        return pl.concat(dfs)

    logging.info("Grouping targets into paired matches and pooled samples for clusters of size 3+")
    sparse_edges = (
        unbinned
        .select(
            "target",
            "coverage",
            samples = pl.col("sample").map_elements(lambda x: [x], return_dtype=pl.List(pl.Utf8)),
            )
        .group_by("target")
        .map_groups(process_groups)
        .group_by(["style", "cluster_size", "samples"])
        .agg(target_ids = pl.col("target").cast(pl.Utf8).sort().str.concat(","))
    )

    return unbinned.with_columns(pl.col("target").cast(pl.Utf8)), sparse_edges

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    MIN_COASSEMBLY_COVERAGE = snakemake.params.min_coassembly_coverage
    MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    PRECLUSTER_SIZE = snakemake.params.precluster_size
    unbinned_path = snakemake.input.unbinned
    distances_path = snakemake.input.distances
    targets_path = snakemake.output.output_targets
    edges_path = snakemake.output.output_edges
    samples = set(snakemake.params.samples)

    unbinned = pl.read_csv(unbinned_path, separator="\t")

    if distances_path:
        sample_distances = (
            pl.scan_csv(distances_path)
            .select("query_name", "match_name", "jaccard")
            .filter(pl.col("jaccard") > 0.01)
            .collect(streaming=True)
        )
        sample_preclusters = get_clusters(
            sample_distances,
            PRECLUSTER_SIZE=PRECLUSTER_SIZE,
            MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
            )
        streaming_pipeline(
            unbinned,
            samples,
            sample_preclusters=sample_preclusters,
            targets_path=targets_path,
            edges_path=edges_path,
            MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
            TAXA_OF_INTEREST=TAXA_OF_INTEREST,
            MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
            CHUNK_SIZE=1000,
            )
    else:
        targets, edges = pipeline(
            unbinned,
            samples,
            MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
            TAXA_OF_INTEREST=TAXA_OF_INTEREST,
            MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
            )
        targets.write_csv(targets_path, separator="\t")
        edges.sort("style", "cluster_size", "samples").write_csv(edges_path, separator="\t")
