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


if __name__ == '__main__':
    unittest.main()
