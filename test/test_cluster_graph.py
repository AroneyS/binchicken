#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.cluster_graph import pipeline

ELUSIVE_EDGES_COLUMNS={
    "samples": str,
    "weight": int,
    "target_ids": str,
    }
READ_SIZE_COLUMNS=["sample", "read_size"]
ELUSIVE_CLUSTERS_COLUMNS=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples", "coassembly"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def test_cluster_graph(self):
        elusive_edges = pl.DataFrame([
            ["sample_2.1,sample_1.1", 2, "0,1"],
            ["sample_1.1,sample_3.1", 1, "2"],
        ], schema=ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["sample_1", 1000],
            ["sample_2", 2000],
            ["sample_3", 3000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["sample_1,sample_2", 2, 2, 2, 3000, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_components(self):
        elusive_edges = pl.DataFrame([
            ["1,2", 1, "some"],
            ["1,3", 2, "some"],
            ["2,3", 3, "some"],
            ["4,5", 4, "some"],
            ["4,6", 5, "some"],
            ["5,6", 6, "some"],
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
            ["5,6", 2, 6, 1, 2000, "4,5,6", "coassembly_0"],
            ["2,3", 2, 3, 1, 2000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_single_bud(self):
        elusive_edges = pl.DataFrame([
            ["1,2", 2, "1,2"],
            ["1,3", 2, "1,3"],
            ["1,4", 2, "1,4"],
            ["2,3", 2, "2,3"],
            ["2,4", 2, "2,4"],
            ["3,4", 2, "3,4"],
            ["4,5", 1, "5"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["5", 1, 0, 0, 1000, "4,5", "coassembly_0"],
            ["4", 1, 0, 0, 1000, "1,2,3,4", "coassembly_1"],
            ["3", 1, 0, 0, 1000, "1,2,3,4", "coassembly_2"],
            ["2", 1, 0, 0, 1000, "1,2,3,4", "coassembly_3"],
            ["1", 1, 0, 0, 1000, "1,2,3,4", "coassembly_4"],
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
            ["1,2", 4, "1,2,3,4"],
            ["3,1", 1, "2"],
            ["3,2", 2, "1,3"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["3", 1, 0, 0, 1000, "2,3", "coassembly_0"],
            ["2", 1, 0, 0, 1000, "1,2", "coassembly_1"],
            ["1", 1, 0, 0, 1000, "1,2", "coassembly_2"],
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
            ["1,2", 3, "1,2,3"],
            ["1,3", 2, "1,3"],
            ["1,4", 2, "1,4"],
            ["2,3", 2, "2,3"],
            ["2,4", 2, "2,4"],
            ["3,4", 3, "1,3,4"],
            ["4,5", 1, "5"],
            ["4,6", 1, "5"],
            ["5,6", 3, "5,6,7"],
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
            ["5,6", 2, 3, 3, 2000, "4,5,6", "coassembly_0"],
            ["3,4", 2, 3, 3, 2000, "1,2,3,4", "coassembly_1"],
            ["1,2", 2, 3, 3, 2000, "1,2,3,4", "coassembly_2"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_choice(self):
        elusive_edges = pl.DataFrame([
            ["1,2", 3, "1,2,3"],
            ["1,3", 2, "1,3"],
            ["2,3", 2, "1,3"],
            ["4,1", 1, "4"],
            ["4,3", 1, "5"],
            ["5,1", 1, "4"],
            ["5,3", 1, "6"],
            ["4,5", 3, "4,5,6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["4,5", 2, 3, 3, 2000, "3,4,5", "coassembly_0"],
            ["1,2", 2, 3, 3, 2000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_irrelevant_targets(self):
        elusive_edges = pl.DataFrame([
            ["1,2", 3, "1,2,3"],
            ["1,3", 2, "1,3"],
            ["2,3", 2, "1,3"],
            ["4,1", 1, "4"],
            ["4,3", 1, "7"],
            ["5,1", 1, "4"],
            ["5,3", 1, "8"],
            ["4,5", 3, "4,5,6"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["4,5", 2, 3, 3, 2000, "3,4,5", "coassembly_0"],
            ["1,2", 2, 3, 3, 2000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_samples_among_many(self):
        elusive_edges = pl.DataFrame([
            ["1,2", 1, "some"],
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
            ["1,2", 2, 1, 1, 2000, "1,2", "coassembly_0"],
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
            ["1,2", 1, "some"],
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
        # 6:               8 9 10 11 12

        elusive_edges = pl.DataFrame([
            ["1,2", 3, "1,2,3"],
            ["1,3", 2, "1,3"],
            ["2,3", 2, "1,3"],
            ["4,1", 1, "4"],
            ["4,3", 1, "5"],
            ["4,5", 2, "6,7"],
            ["4,6", 2, "8,9"],
            ["5,6", 3, "10,11,12"],
            ["1,2,3", 2, "1,3"],
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
            ["4,5,6", 3, 7, 7, 3000, "4,5,6", "coassembly_0"],
            ["1,2,3", 3, 7, 3, 3000, "1,2,3", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            MIN_COASSEMBLY_SAMPLES=3,
            MAX_COASSEMBLY_SAMPLES=3,
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
            ["1,2", 2, "3,4"],
            ["1,3", 2, "2,4"],
            ["1,4", 3, "2,3,4"],
            ["2,3", 2, "1,4"],
            ["2,4", 3, "1,3,4"],
            ["3,4", 3, "1,2,4"],
            # triplets of 1,2,3,4
            ["1,2,3", 1, "4"],
            ["1,2,4", 2, "3,4"],
            ["1,3,4", 2, "2,4"],
            ["2,3,4", 2, "1,4"],
            # quads of 1,2,3,4
            ["1,2,3,4", 1, "4"],
            # pairs of 5,6,7,8
            ["5,6", 2, "7,8"],
            ["5,7", 2, "6,8"],
            ["5,8", 3, "8,9,10"],
            ["6,7", 2, "5,8"],
            ["6,8", 1, "8"],
            ["7,8", 1, "8"],
            # triplets of 5,6,7,8
            ["5,6,7", 1, "8"],
            ["5,6,8", 1, "8"],
            ["5,7,8", 1, "8"],
            ["6,7,8", 1, "8"],
            # quads of 5,6,7,8
            ["5,6,7,8", 1, "8"],
            # joint pairs
            ["2,5", 1, "1"],
            ["3,5", 1, "1"],
            ["4,5", 1, "1"],
            # joint triplets
            ["2,3,5", 1, "1"],
            ["2,4,5", 1, "1"],
            ["3,4,5", 1, "1"],
            # joint quads
            ["2,3,4,5", 1, "1"],
        ], schema = ELUSIVE_EDGES_COLUMNS)
        read_size = pl.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
            ["7", 1000],
            ["8", 1000],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["5,6,7,8", 4, 9, 7, 4000, "5,6,7,8", "coassembly_0"],
            ["1,2,3,4", 4, 10, 3, 4000, "1,2,3,4", "coassembly_1"],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=4,
            MIN_COASSEMBLY_SAMPLES=4,
            MAX_COASSEMBLY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

if __name__ == '__main__':
    unittest.main()
