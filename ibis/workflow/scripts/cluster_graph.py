########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import polars as pl
import itertools
import os
import logging

OUTPUT_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": int,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

def join_list_subsets(df1, df2, list_colname, value_colname, output_colname):
    """
    Join two DataFrames/LazyFrames by strict subset of list column
    """
    df2 = (
        df2.select(
            pl.col(list_colname),
            pl.col(value_colname),
            )
        .with_row_count()
    )

    output = (
        df1
        .select(
            pl.col(list_colname),
            original_list = pl.col(list_colname),
            )
        .with_row_count()
        .explode(list_colname)
        .join(
            df2
            .explode(list_colname),
            on=list_colname,
        )
        .groupby('row_nr', 'row_nr_right')
        .agg(list_colname, pl.first("original_list"))
        .join(
            df2,
            left_on='row_nr_right', right_on='row_nr'
        )
        .filter(pl.col(list_colname).hash() == pl.col(list_colname + '_right').hash())
        .groupby("original_list")
        .agg(pl.col(value_colname).flatten())
        .select(
            pl.col(value_colname).list.unique().alias(output_colname),
            join_col = pl.col("original_list").list.sort().hash(),
            )
    )

    return (
        df1
        .with_columns(join_col = pl.col(list_colname).list.sort().hash())
        .join(output, on="join_col", how="left")
        .drop("join_col")
        .with_columns(pl.col(output_colname).fill_null([]))
        )

def accumulate_clusters(x):
    clustered_samples = set()
    choices = []
    for cluster in x:
        if all([x not in clustered_samples for x in cluster]):
            choices.append(True)
            clustered_samples.update(cluster)
        else:
            choices.append(False)

    return pl.Series(choices, dtype=pl.Boolean)

def find_recover_candidates(df, samples_df, MAX_RECOVERY_SAMPLES=20):
    samples_df = samples_df.explode("target_ids")
    df = df.with_columns(join_col = pl.col("samples").hash())

    output = (
        df
        .select("join_col", "target_ids")
        .explode("target_ids")
        .join(samples_df, on="target_ids", how="left")
        .groupby("join_col", "recover_candidates")
        .agg(pl.count())
        .sort("count", descending=True)
        .groupby("join_col")
        .head(MAX_RECOVERY_SAMPLES)
        .groupby("join_col")
        .agg("recover_candidates")
    )

    return (
        df
        .join(output, on="join_col", how="left")
        .drop("join_col")
    )

def pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20):

    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    if len(elusive_edges) == 0:
        logging.warning("No elusive edges found")
        return pl.DataFrame(schema=OUTPUT_COLUMNS)

    is_pooled = any(elusive_edges["style"] == "pool")

    with pl.StringCache():
        elusive_edges = (
            elusive_edges
            .lazy()
            .with_columns(
                pl.col("samples")
                    .str.split(",")
                    .list.eval(pl.element().str.replace(r"\.1$", ""))
                    .cast(pl.List(pl.Categorical)),
                pl.col("target_ids")
                    .str.split(",")
                    .cast(pl.List(pl.UInt32)),
                )
            .with_columns(length = pl.col("samples").list.lengths())
        )

        read_size = (
            read_size
            .lazy()
            .with_columns(
                pl.col("sample").cast(pl.Categorical),
                )
        )

        if MAX_COASSEMBLY_SAMPLES == 1:
            logging.info("Skipping clustering, using single-sample clusters")
            clusters = [
                elusive_edges
                .explode("samples")
                .groupby("samples")
                .agg(pl.col("target_ids").flatten())
                .with_columns(
                    pl.col("samples").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)).cast(pl.List(pl.Categorical)),
                    pl.col("target_ids").list.sort().list.unique(),
                    sample = pl.col("samples").apply(lambda x: [x], return_dtype=pl.List(pl.Utf8)).cast(pl.List(pl.Categorical)),
                )
            ]
        else:
            logging.info("Forming candidate sample clusters")
            clusters = [
                # Pair clusters
                elusive_edges
                .filter(pl.col("style") == "match")
                .filter(pl.col("length") == 2)
                .select("samples", "target_ids", sample = pl.col("samples"))
            ]

            if is_pooled:
                    clusters.append(
                        # 3+ way clusters
                        elusive_edges
                        .filter(pl.col("style") == "pool")
                        # Prevent combinatorial explosion (also, large clusters are less useful for distinguishing between clusters)
                        .filter(pl.col("samples").list.lengths() < 100)
                        .with_columns(
                            samples_combinations = pl.struct(["length", "cluster_size"]).apply(
                                lambda x: [i for i in itertools.combinations(range(x["length"]), x["cluster_size"])],
                                return_dtype=pl.List(pl.List(pl.Int64)),
                                )
                            )
                        .explode("samples_combinations")
                        .with_columns(
                            pl.col("samples")
                            .list.take(pl.col("samples_combinations"))
                            )
                        .select("samples", "target_ids")
                        .groupby("samples")
                        .agg(pl.col("target_ids").flatten())
                        .pipe(
                            join_list_subsets,
                            df2=elusive_edges,
                            list_colname="samples",
                            value_colname="target_ids",
                            output_colname="extra_targets"
                            )
                        .select(
                            "samples",
                            pl.concat_list("target_ids", "extra_targets").list.unique(),
                            sample = pl.col("samples"),
                            )
                    )

        sample_targets = (
            elusive_edges
            .select("target_ids", recover_candidates = pl.col("samples"))
            .explode("recover_candidates")
            .explode("target_ids")
            .unique()
            .groupby("recover_candidates")
            .agg("target_ids")
        )

        logging.info("Filtering clusters (each sample restricted to only once per cluster size)")
        clusters = (
            pl.concat(clusters)
            .with_columns(
                "samples", "target_ids",
                total_targets = pl.col("target_ids").list.lengths(),
                length = pl.col("samples").list.lengths()
                )
            .filter(pl.col("length") >= MIN_COASSEMBLY_SAMPLES)
            .explode("sample")
            .join(read_size, on="sample", how="inner")
            .groupby("samples")
            .agg(
                pl.first("length"),
                pl.first("target_ids"),
                pl.first("total_targets"),
                total_size = pl.sum("read_size"),
                )
            .filter(
                pl.when(MAX_COASSEMBLY_SIZE is not None)
                .then(pl.col("total_size") <= MAX_COASSEMBLY_SIZE)
                .otherwise(True)
                )
            .sort("total_targets", descending=True)
            .with_columns(
                unique_samples = 
                    pl.col("samples")
                    .map(accumulate_clusters, return_dtype=pl.Boolean),
                )
            .filter(pl.col("unique_samples"))
            .pipe(
                find_recover_candidates,
                sample_targets,
                MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES,
                )
            .with_columns(
                pl.col("samples").cast(pl.List(pl.Utf8)).list.sort().list.join(","),
                recover_samples = pl.concat_list(pl.col("samples"), "recover_candidates")
                .list.unique(maintain_order=True)
                .list.head(MAX_RECOVERY_SAMPLES)
                .cast(pl.List(pl.Utf8))
                .list.sort()
                .list.join(",")
                )
            .sort("total_targets", "samples", descending=True)
            .with_row_count("coassembly")
            .select(
                "samples", "length", "total_targets", "total_size", "recover_samples",
                coassembly = pl.lit("coassembly_") + pl.col("coassembly").cast(pl.Utf8)
                )
            .collect(streaming=True)
        )

    logging.info("Done")

    return clusters

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    MAX_COASSEMBLY_SIZE = snakemake.params.max_coassembly_size * 10**9 if snakemake.params.max_coassembly_size else None
    MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
    MIN_COASSEMBLY_SAMPLES = snakemake.params.num_coassembly_samples
    MAX_RECOVERY_SAMPLES = snakemake.params.max_recovery_samples
    elusive_edges_path = snakemake.input.elusive_edges
    read_size_path = snakemake.input.read_size
    elusive_clusters_path = snakemake.output.elusive_clusters

    elusive_edges = pl.read_csv(elusive_edges_path, separator="\t", dtypes={"target_ids": str})
    read_size = pl.read_csv(read_size_path, has_header=False, new_columns=["sample", "read_size"])

    clusters = pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES
        )
    clusters.write_csv(elusive_clusters_path, separator="\t")
