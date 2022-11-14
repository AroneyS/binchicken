#!/usr/bin/env python3

import unittest
import pandas as pd
from cockatoo.workflow.scripts.cluster_graph import pipeline

ELUSIVE_EDGES_COLUMNS=["taxa_group", "weight", "target_ids", "sample1", "sample2"]
READ_SIZE_COLUMNS=["sample", "read_size"]
ELUSIVE_CLUSTERS_COLUMNS=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples", "coassembly"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def test_cluster_graph(self):
        elusive_edges = pd.DataFrame([
            ["Root", 2, "0,1", "sample_2.1", "sample_1.1"],
            ["Root", 1, "2", "sample_1.1", "sample_3.1"],
        ], columns=ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["sample_1", 1000],
            ["sample_2", 2000],
            ["sample_3", 3000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["sample_1,sample_2", 2, 2, 2, 3000, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_components(self):
        elusive_edges = pd.DataFrame([
            ["Root", 1, "some", "1", "2"],
            ["Root", 2, "some", "1", "3"],
            ["Root", 3, "some", "2", "3"],
            ["Root", 4, "some", "4", "5"],
            ["Root", 5, "some", "4", "6"],
            ["Root", 6, "some", "5", "6"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["5,6", 2, 6, 1, 2000, "4,5,6", "coassembly_0"],
            ["2,3", 2, 3, 1, 2000, "1,2,3", "coassembly_1"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_single_bud(self):
        elusive_edges = pd.DataFrame([
            ["Root", 2, "1,2", "1", "2"],
            ["Root", 2, "1,3", "1", "3"],
            ["Root", 2, "1,4", "1", "4"],
            ["Root", 2, "2,3", "2", "3"],
            ["Root", 2, "2,4", "2", "4"],
            ["Root", 2, "3,4", "3", "4"],
            ["Root", 1, "5", "4", "5"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["5", 1, 0, 0, 1000, "4,5", "coassembly_0"],
            ["4", 1, 0, 0, 1000, "1,2,3,4", "coassembly_1"],
            ["3", 1, 0, 0, 1000, "1,2,3,4", "coassembly_2"],
            ["2", 1, 0, 0, 1000, "1,2,3,4", "coassembly_3"],
            ["1", 1, 0, 0, 1000, "1,2,3,4", "coassembly_4"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_COASSEMBLY_SAMPLES=1,
            MIN_COASSEMBLY_SAMPLES=1,
            MAX_RECOVERY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_single_bud_choice(self):
        elusive_edges = pd.DataFrame([
            ["Root", 4, "1,2,3,4", "1", "2"],
            ["Root", 1, "2", "3", "1"],
            ["Root", 2, "1,3", "3", "2"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["3", 1, 0, 0, 1000, "2,3", "coassembly_0"],
            ["2", 1, 0, 0, 1000, "1,2", "coassembly_1"],
            ["1", 1, 0, 0, 1000, "1,2", "coassembly_2"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_COASSEMBLY_SAMPLES=1,
            MIN_COASSEMBLY_SAMPLES=1,
            MAX_RECOVERY_SAMPLES=2,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud(self):
        elusive_edges = pd.DataFrame([
            ["Root", 3, "1,2,3", "1", "2"],
            ["Root", 2, "1,3", "1", "3"],
            ["Root", 2, "1,4", "1", "4"],
            ["Root", 2, "2,3", "2", "3"],
            ["Root", 2, "2,4", "2", "4"],
            ["Root", 3, "1,3,4", "3", "4"],
            ["Root", 1, "5", "4", "5"],
            ["Root", 1, "5", "4", "6"],
            ["Root", 3, "5,6,7", "5", "6"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["5,6", 2, 3, 3, 2000, "4,5,6", "coassembly_0"],
            ["3,4", 2, 3, 3, 2000, "1,2,3,4", "coassembly_1"],
            ["1,2", 2, 3, 3, 2000, "1,2,3,4", "coassembly_2"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=4,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_choice(self):
        elusive_edges = pd.DataFrame([
            ["Root", 3, "1,2,3", "1", "2"],
            ["Root", 2, "1,3", "1", "3"],
            ["Root", 2, "1,3", "2", "3"],
            ["Root", 1, "4", "4", "1"],
            ["Root", 1, "5", "4", "3"],
            ["Root", 1, "4", "5", "1"],
            ["Root", 1, "6", "5", "3"],
            ["Root", 3, "4,5,6", "4", "5"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["4,5", 2, 3, 3, 2000, "3,4,5", "coassembly_0"],
            ["1,2", 2, 3, 3, 2000, "1,2,3", "coassembly_1"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_double_bud_irrelevant_targets(self):
        elusive_edges = pd.DataFrame([
            ["Root", 3, "1,2,3", "1", "2"],
            ["Root", 2, "1,3", "1", "3"],
            ["Root", 2, "1,3", "2", "3"],
            ["Root", 1, "4", "4", "1"],
            ["Root", 1, "7", "4", "3"],
            ["Root", 1, "4", "5", "1"],
            ["Root", 1, "8", "5", "3"],
            ["Root", 3, "4,5,6", "4", "5"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["4,5", 2, 3, 3, 2000, "1,4,5", "coassembly_0"],
            ["1,2", 2, 3, 3, 2000, "1,2,3", "coassembly_1"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(
            elusive_edges,
            read_size,
            MAX_RECOVERY_SAMPLES=3,
            )
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_two_samples_among_many(self):
        elusive_edges = pd.DataFrame([
            ["Root", 1, "some", "1", "2"],
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
            ["1,2", 2, 1, 1, 2000, "1,2", "coassembly_0"],
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_cluster_no_edges(self):
        elusive_edges = pd.DataFrame([
        ], columns = ELUSIVE_EDGES_COLUMNS)
        read_size = pd.DataFrame([
            ["1", 1000],
            ["2", 1000],
            ["3", 1000],
            ["4", 1000],
            ["5", 1000],
            ["6", 1000],
        ], columns=READ_SIZE_COLUMNS)

        expected = pd.DataFrame([
        ], columns=ELUSIVE_CLUSTERS_COLUMNS).astype(object)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)

if __name__ == '__main__':
    unittest.main()
