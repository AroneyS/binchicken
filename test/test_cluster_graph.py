#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal, assert_series_equal
from ibis.workflow.scripts.cluster_graph import pipeline, join_list_subsets, accumulate_clusters, find_recover_candidates

ELUSIVE_EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
    }
READ_SIZE_COLUMNS=["sample", "read_size"]
ELUSIVE_CLUSTERS_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": int,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }

CAT_CLUSTERS_COLUMNS={
    "samples_hash": pl.UInt64,
    "samples": pl.List(pl.Utf8),
    "length": int,
    "target_ids": pl.List(pl.UInt32),
    "total_targets": int,
    "total_size": int,
    }
SAMPLE_TARGETS_COLUMNS={
    "recover_candidates": pl.Utf8,
    "target_ids": pl.List(pl.UInt32),
    }
CAT_RECOVERY_COLUMNS={
    "samples_hash": pl.UInt64,
    "samples": pl.List(pl.Utf8),
    "length": int,
    "target_ids": pl.List(pl.UInt32),
    "total_targets": int,
    "total_size": int,
    "recover_candidates": pl.List(pl.Utf8),
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def assertSeriesEqual(self, a, b):
        assert_series_equal(a, b, check_dtype=False)

    def test_cluster_graph(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "sample_2.1,sample_1.1", "0,1,2"],
            ["match", 2, "sample_1.1,sample_3.1", "1,2"],
        ], schema=ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["sample_1", 1000],
            ["sample_2", 2000],
            ["sample_3", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["sample_1,sample_2", 2, 3, 3000, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_components(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1"],
            ["match", 2, "1,3", "1,2"],
            ["match", 2, "2,3", "1,2,3"],
            ["match", 2, "4,5", "4,5,6,7"],
            ["match", 2, "4,6", "4,5,6,7,8"],
            ["match", 2, "5,6", "4,5,6,7,8,9"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["5,6", 2, 6, 2000, "4,5,6", "coassembly_0"],
            ["2,3", 2, 3, 2000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_single_bud(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "1,4", "1,4"],
            ["match", 2, "2,3", "2,3"],
            ["match", 2, "2,4", "2,4"],
            ["match", 2, "3,4", "3,4"],
            ["match", 2, "4,5", "5"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["4", 1, 5, 4000, "1,2,3,4", "coassembly_0"],
            ["1", 1, 4, 1000, "1,2,3,4", "coassembly_1"],
            ["2", 1, 4, 2000, "1,2,3,4", "coassembly_2"],
            ["3", 1, 4, 3000, "1,2,3,4", "coassembly_3"],
            ["5", 1, 1, 5000, "4,5", "coassembly_4"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_COASSEMBLY_SAMPLES=1,
            MIN_COASSEMBLY_SAMPLES=1,
            MAX_RECOVERY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_single_bud_choice(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3,4"],
            ["match", 2, "3,1", "5"],
            ["match", 2, "3,2", "6,7"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["2", 1, 6, 1000, "1,2", "coassembly_0"],
            ["1", 1, 5, 1000, "1,2", "coassembly_1"],
            ["3", 1, 3, 1000, "2,3", "coassembly_2"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_COASSEMBLY_SAMPLES=1,
            MIN_COASSEMBLY_SAMPLES=1,
            MAX_RECOVERY_SAMPLES=2,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "1,4", "1,4"],
            ["match", 2, "2,3", "2,3"],
            ["match", 2, "2,4", "2,4"],
            ["match", 2, "3,4", "1,3,4"],
            ["match", 2, "4,5", "5"],
            ["match", 2, "4,6", "5"],
            ["match", 2, "5,6", "5,6,7"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
            ["6", 6000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2", 2, 3, 3000, "1,2,3,4", "coassembly_0"],
            ["3,4", 2, 3, 7000, "1,2,3,4", "coassembly_1"],
            ["5,6", 2, 3, 11000, "4,5,6", "coassembly_2"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_choice(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "2,3", "1,3"],
            ["match", 2, "4,1", "4"],
            ["match", 2, "4,3", "5"],
            ["match", 2, "5,1", "4"],
            ["match", 2, "5,3", "6"],
            ["match", 2, "4,5", "4,5,6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2", 2, 3, 3000, "1,2,3", "coassembly_0"],
            ["4,5", 2, 3, 9000, "3,4,5", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_irrelevant_targets(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "2,3", "1,3"],
            ["match", 2, "4,1", "4"],
            ["match", 2, "4,3", "7"],
            ["match", 2, "5,1", "4"],
            ["match", 2, "5,3", "8"],
            ["match", 2, "4,5", "4,5,6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2", 2, 3, 3000, "1,2,3", "coassembly_0"],
            ["4,5", 2, 3, 9000, "1,4,5", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_samples_among_many(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "9,10"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2", 2, 2, 2000, "1,2", "coassembly_0"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_no_edges(self):
        elusive_edges = pl.DataFrame([
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_only_large_clusters(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "9,10"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 10000],
            ["2", 10000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size, MAX_COASSEMBLY_SIZE=2000)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_three_samples(self):
        # 1: 1 2 3 4
        # 2: 1 2 3
        # 3: 1   3   5
        # 4:       4 5 6 7 8 9
        # 5:           6 7     10 11 12
        # 6:           6   8 9 10 11 12

        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "2,3", "1,3"],
            ["match", 2, "4,1", "4"],
            ["match", 2, "4,3", "5"],
            ["match", 2, "4,5", "6,7"],
            ["match", 2, "4,6", "8,9"],
            ["match", 2, "5,6", "10,11,12"],
            ["pool", 3, "1,2,3", "1,3"],
            ["pool", 3, "4,5,6", "6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
            ["6", 6000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2,3", 3, 2, 6000, "1,2,3", "coassembly_0"],
            ["4,5,6", 3, 1, 15000, "4,5,6", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            MIN_COASSEMBLY_SAMPLES=3,
            MAX_COASSEMBLY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_three_samples_min_cluster(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1,2,3"],
            ["match", 2, "1,3", "1,3"],
            ["match", 2, "2,3", "1,3"],
            ["match", 2, "4,1", "4"],
            ["match", 2, "4,3", "5"],
            ["match", 2, "4,5", "6,7"],
            ["match", 2, "4,6", "8,9"],
            ["match", 2, "5,6", "10,11,12"],
            ["pool", 3, "1,2,3", "1,3"],
            ["pool", 3, "4,5,6", "6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
            ["6", 6000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2,3", 3, 2, 6000, "1,2,3", "coassembly_0"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            MIN_COASSEMBLY_SAMPLES=3,
            MAX_COASSEMBLY_SAMPLES=3,
            MIN_CLUSTER_TARGETS=2,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_three_samples_max_sample(self):
        elusive_edges = pl.DataFrame([
            ["pool", 3, "1,2,3", "1,3"],
            ["pool", 3, "1,2,3,4,5,6", "2"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
            ["6", 6000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,2,3", 3, 3, 6000, "1,2,3,4,5,6", "coassembly_0"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=6,
            MIN_COASSEMBLY_SAMPLES=3,
            MAX_COASSEMBLY_SAMPLES=3,
            MAX_SAMPLES_COMBINATIONS=5,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_four_samples(self):
        # 1:   2 3 4
        # 2: 1   3 4
        # 3: 1 2   4
        # 4: 1 2 3 4

        # 5: 1         6 7 8 9 10
        # 6:         5   7 8
        # 7:         5 6   8
        # 8:               8 9 10

        elusive_edges = pl.DataFrame([
            # pairs of 1,2,3,4
            ["match", 2, "1,2", "3,4"],
            ["match", 2, "1,3", "2,4"],
            ["match", 2, "1,4", "2,3,4"],
            ["match", 2, "2,3", "1,4"],
            ["match", 2, "2,4", "1,3,4"],
            ["match", 2, "3,4", "1,2,4"],
            # pairs of 5,6,7,8
            ["match", 2, "5,6", "7,8"],
            ["match", 2, "5,7", "6,8"],
            ["match", 2, "5,8", "8,9,10"],
            ["match", 2, "6,7", "5,8"],
            ["match", 2, "6,8", "8"],
            ["match", 2, "7,8", "8"],
            # joint pairs
            ["match", 2, "2,5", "1"],
            ["match", 2, "3,5", "1"],
            ["match", 2, "4,5", "1"],
            # triplets
            ["pool", 3, "2,3,4,5", "1"],
            ["pool", 3, "1,3,4", "2"],
            ["pool", 3, "1,2,4", "3"],
            ["pool", 3, "1,2,3,4", "4"],
            ["pool", 3, "5,6,7,8", "8"],
            # quads
            ["pool", 4, "2,3,4,5", "1"],
            ["pool", 4, "1,2,3,4", "4"],
            ["pool", 4, "5,6,7,8", "8"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 2000],
            ["3", 3000],
            ["4", 4000],
            ["5", 5000],
            ["6", 6000],
            ["7", 7000],
            ["8", 8000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["1,4", 2, 3, 5000, "1,2,3,4", "coassembly_0"],
            ["5,8", 2, 3, 13000, "5,6,7,8", "coassembly_1"],
            ["2,3", 2, 2, 5000, "2,3,4,5", "coassembly_2"],
            ["1,2,4", 3, 2, 7000, "1,2,3,4", "coassembly_3"],
            ["6,7", 2, 2, 13000, "5,6,7,8", "coassembly_4"],
            ["1,2,3,4", 4, 1, 10000, "1,2,3,4", "coassembly_5"],
            ["5,6,7", 3, 1, 18000, "5,6,7,8", "coassembly_6"],
            ["5,6,7,8", 4, 1, 26000, "5,6,7,8", "coassembly_7"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=4,
            MIN_COASSEMBLY_SAMPLES=2,
            MAX_COASSEMBLY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_exclude_coassemblies(self):
        elusive_edges = pl.DataFrame([
            ["match", 2, "1,2", "1"],
            ["match", 2, "1,3", "1,2"],
            ["match", 2, "2,3", "1,2,3"],
            ["match", 2, "4,5", "4,5,6,7"],
            ["match", 2, "4,6", "4,5,6,7,8"],
            ["match", 2, "5,6", "4,5,6,7,8,9"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["4,6", 2, 5, 2000, "4,5,6", "coassembly_0"],
            ["1,3", 2, 2, 2000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size, EXCLUDE_COASSEMBLIES=["2,3", "5,6"])
        self.assertDataFrameEqual(expected, observed)

    def test_join_list_subsets(self):
        with pl.StringCache():
            df1 = (
                pl.DataFrame([
                        [["a", "b"], ["1"], 2],
                        [["b", "c", "d"], ["2"], 3],
                        [["a", "b", "c", "d"], ["3"], 4],
                        [["d", "e", "f"], ["4"], 3],
                        [["b", "c", "d", "g"], ["5"], 4],
                    ], schema=["samples", "target_ids", "length"])
                .with_columns(
                    pl.col("samples").cast(pl.List(pl.Categorical)),
                    pl.col("length").cast(pl.UInt32),
                    )
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
            )

            df2 = (
                pl.DataFrame([
                        [["a", "b"], ["2"], 2],
                        [["b", "c"], ["2"], 2],
                        [["c", "d"], ["3"], 2],
                        [["b", "c", "d"], ["4"], 3],
                        [["b", "c", "d", "g"], ["5"], 3],
                        [["b", "c", "d"], ["6"], 2],
                    ], schema=["samples", "target_ids", "cluster_size"])
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
            )

            expected = (
                pl.DataFrame([
                        [["a", "b"], ["1"], 2, ["2"]],
                        [["b", "c", "d"], ["2"], 3, ["4", "5"]],
                        [["a", "b", "c", "d"], ["3"], 4, []],
                        [["d", "e", "f"], ["4"], 3, []],
                        [["b", "c", "d", "g"], ["5"], 4, []],
                    ], schema=["samples", "target_ids", "length", "extra_targets"])
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
                .select("samples", "target_ids", "length", "samples_hash", "extra_targets")
            )

            observed = (
                join_list_subsets(df1, df2)
                .with_columns(pl.col("extra_targets").list.sort())
            )
            self.assertDataFrameEqual(expected, observed)

    def test_join_list_subsets_lazy(self):
        with pl.StringCache():
            df1 = (
                pl.DataFrame([
                        [["a", "b"], ["1"], 2],
                        [["b", "c", "d"], ["2"], 3],
                        [["a", "b", "c", "d"], ["3"], 4],
                        [["d", "e", "f"], ["4"], 3],
                        [["b", "c", "d", "g"], ["5"], 4],
                    ], schema=["samples", "target_ids", "length"])
                .lazy()
                .with_columns(
                    pl.col("samples").cast(pl.List(pl.Categorical)),
                    pl.col("length").cast(pl.UInt32),
                    )
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
            )

            df2 = (
                pl.DataFrame([
                        [["a", "b"], ["2"], 2],
                        [["b", "c"], ["2"], 2],
                        [["c", "d"], ["3"], 2],
                        [["b", "c", "d"], ["4"], 3],
                        [["b", "c", "d", "g"], ["5"], 3],
                        [["b", "c", "d"], ["6"], 2],
                    ], schema=["samples", "target_ids", "cluster_size"])
                .lazy()
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
            )

            expected = (
                pl.DataFrame([
                        [["a", "b"], ["1"], 2, ["2"]],
                        [["b", "c", "d"], ["2"], 3, ["4", "5"]],
                        [["a", "b", "c", "d"], ["3"], 4, []],
                        [["d", "e", "f"], ["4"], 3, []],
                        [["b", "c", "d", "g"], ["5"], 4, []],
                    ], schema=["samples", "target_ids", "length", "extra_targets"])
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
                .with_columns(samples_hash = pl.col("samples").list.sort().hash())
                .select("samples", "target_ids", "length", "samples_hash", "extra_targets")
            )

            observed = (
                join_list_subsets(df1, df2)
                .with_columns(pl.col("extra_targets").list.sort())
                .collect()
            )
            self.assertDataFrameEqual(expected, observed)

    def test_accumulate_clusters(self):
        input = [[1,2,3], [1,2,4], [4,5,6], [4,5,7], [7,8,9], [7,8,10], [10,11,12], [10,11]]

        expected = pl.Series([True, False, True, False, True, False, True, True], dtype=pl.Boolean)
        observed = accumulate_clusters(input)
        self.assertSeriesEqual(observed, expected)

    def test_accumulate_clusters_diff_sample_sizes(self):
        input = [[1,2,3], [1,2], [1,2,4], [2,3], [4,5,6], [3,4], [4,5,7], [4,5], [7,8,9], [7,8,10], [10,11,12], [10,11], [1,2,3,4], [1,2,3,4,5], [1,2,3,4,5,6]]

        expected = pl.Series([True, True, False, False, True, True, False, False, True, False, True, True, True, True, True], dtype=pl.Boolean)
        observed = accumulate_clusters(input)
        self.assertSeriesEqual(observed, expected)

    def test_accumulate_clusters_empty_input(self):
        input = []

        expected = pl.Series([], dtype=pl.Boolean)
        observed = accumulate_clusters(input)
        self.assertSeriesEqual(observed, expected)

    def test_find_recover_candidates(self):
        with pl.StringCache():
            sample_targets = (
                pl.DataFrame([
                        ["1",  [1,2,3]],
                        ["2",  [1,2,3]],
                        ["3",  [1,2,3]],
                        ["4",  [1,2,3]],
                        ["5",  [4,5,6]],
                        ["6",  [4,5,6]],
                        ["7",  [4,5,6]],
                        ["8",  [4,5,6]],
                        ["9",  [1,2,3]],
                        ["10", [1,2]],
                        ["11", [1]],
                        ["12", [4,5,6]],
                        ["13", [4,5]],
                        ["14", [4]],
                    ],
                    schema=SAMPLE_TARGETS_COLUMNS
                    )
                .with_columns(pl.col("recover_candidates").cast(pl.Categorical))
            )

            clusters = (
                pl.DataFrame([
                        [12341234, ["1","2","3","4"], 4, [1,2,3], 3, 4000],
                        [56785678, ["5","6","7","8"], 4, [4,5,6], 3, 4000],
                    ],
                    schema=CAT_CLUSTERS_COLUMNS
                    )
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
            )

            expected = (
                pl.DataFrame([
                        [12341234, ["1","2","3","4"], 4, [1,2,3], 3, 4000, ["1","2","3","4","9","10"]],
                        [56785678, ["5","6","7","8"], 4, [4,5,6], 3, 4000, ["5","6","7","8","12","13"]],
                    ],
                    schema=CAT_RECOVERY_COLUMNS
                    )
                .with_columns(
                    pl.col("samples").cast(pl.List(pl.Categorical)),
                    pl.col("recover_candidates").cast(pl.List(pl.Categorical)),
                    )
            )

            observed = (
                find_recover_candidates(
                    clusters,
                    sample_targets,
                    MAX_RECOVERY_SAMPLES=6,
                    )
                .with_columns(pl.col("recover_candidates").list.sort())
            )
            self.assertDataFrameEqual(expected, observed)

    def test_find_recover_candidates_lazy(self):
        with pl.StringCache():
            sample_targets = (
                pl.DataFrame([
                        ["1",  [1,2,3]],
                        ["2",  [1,2,3]],
                        ["3",  [1,2,3]],
                        ["4",  [1,2,3]],
                        ["5",  [4,5,6]],
                        ["6",  [4,5,6]],
                        ["7",  [4,5,6]],
                        ["8",  [4,5,6]],
                        ["9",  [1,2,3]],
                        ["10", [1,2]],
                        ["11", [1]],
                        ["12", [4,5,6]],
                        ["13", [4,5]],
                        ["14", [4]],
                    ],
                    schema=SAMPLE_TARGETS_COLUMNS
                    )
                .with_columns(pl.col("recover_candidates").cast(pl.Categorical))
                .lazy()
            )

            clusters = (
                pl.DataFrame([
                        [12341234, ["1","2","3","4"], 4, [1,2,3], 3, 4000],
                        [56785678, ["5","6","7","8"], 4, [4,5,6], 3, 4000],
                    ],
                    schema=CAT_CLUSTERS_COLUMNS
                    )
                .with_columns(pl.col("samples").cast(pl.List(pl.Categorical)))
                .lazy()
            )

            expected = (
                pl.DataFrame([
                        [12341234, ["1","2","3","4"], 4, [1,2,3], 3, 4000, ["1","2","3","4","9","10"]],
                        [56785678, ["5","6","7","8"], 4, [4,5,6], 3, 4000, ["5","6","7","8","12","13"]],
                    ],
                    schema=CAT_RECOVERY_COLUMNS
                    )
                .with_columns(
                    pl.col("samples").cast(pl.List(pl.Categorical)),
                    pl.col("recover_candidates").cast(pl.List(pl.Categorical)),
                    )
            )

            observed = (
                find_recover_candidates(
                    clusters,
                    sample_targets,
                    MAX_RECOVERY_SAMPLES=6,
                    )
                .with_columns(pl.col("recover_candidates").list.sort())
                .collect(streaming=True)
            )
            self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
