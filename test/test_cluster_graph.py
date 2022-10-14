#!/usr/bin/env python3

import unittest
import pandas as pd
from cockatoo.workflow.scripts.cluster_graph import pipeline

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def test_cluster_graph(self):
        elusive_edges = pd.DataFrame({
            "taxa_group": ["Root","Root"],
            "weight": [2, 1],
            "target_ids": ["0,1", "2"],
            "sample1": ["sample_2.1", "sample_1.1"],
            "sample2": ["sample_1.1", "sample_3.1"],
        })
        read_size = pd.DataFrame({
            "sample": ["sample_1", "sample_2", "sample_3"],
            "read_size": [1000, 2000, 3000],
        })

        expected = pd.DataFrame({
            "samples": ["sample_1,sample_2"],
            "length": [2],
            "total_weight": [2],
            "total_targets": [2],
            "total_size": [3000],
            "recover_samples": ["sample_1,sample_2"],
            "coassembly": ["coassembly_0"],
        }).astype(object)
        observed = pipeline(elusive_edges, read_size)
        self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
