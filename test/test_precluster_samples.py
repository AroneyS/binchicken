#!/usr/bin/env python3

import unittest
import numpy as np
from binchicken.workflow.scripts.precluster_samples import processing

class Tests(unittest.TestCase):
    def assertListListEqual(self, a, b):
        a_sets = [sorted(x) for x in a]
        b_sets = [sorted(x) for x in b]
        self.assertEqual(sorted(a_sets), sorted(b_sets))

    def test_precluster_samples(self):
        distances = 1 - np.array([
            [1.0,  1.0,  0.0, 0.25],
            [1.0,  1.0,  0.0, 0.25],
            [0.0,  0.0,  1.0, 0.5],
            [0.25, 0.25, 0.5, 1.0],
        ])
        samples = ["sample_1", "sample_2", "sample_3", "sample_4"]

        expected = [
            ["sample_1", "sample_2"],
            ["sample_3", "sample_4"],
        ]

        observed = processing(distances, samples, MAX_CLUSTER_SIZE=2)
        self.assertListListEqual(expected, observed)

    def test_precluster_samples_link(self):
        # 1+2 cluster tightly
        # 3 clusters with 1+2
        # 4+5 cluster more loose than 3 with 1+2
        distances = 1 - np.array([
            [1.0, 1.0, 0.9, 0.0, 0.0],
            [1.0, 1.0, 0.9, 0.0, 0.0],
            [0.9, 0.9, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.7],
            [0.0, 0.0, 0.0, 0.7, 1.0],
        ])
        samples = ["sample_1", "sample_2", "sample_3", "sample_4", "sample_5"]

        expected = [
            ["sample_1", "sample_2"],
            ["sample_3"],
            ["sample_4", "sample_5"],
        ]

        observed = processing(distances, samples, MAX_CLUSTER_SIZE=2)
        self.assertListListEqual(expected, observed)

    def test_precluster_samples_singletons(self):
        # 1+2 cluster tightly
        # 3 clusters with 1+2
        # 4+5 cluster more loose than 3 with 1+2
        # 6 is a singleton
        distances = 1 - np.array([
            [1.0, 1.0, 0.9, 0.0, 0.0, 0.0],
            [1.0, 1.0, 0.9, 0.0, 0.0, 0.0],
            [0.9, 0.9, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.7, 0.0],
            [0.0, 0.0, 0.0, 0.7, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ])
        samples = ["sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6"]

        expected = [
            ["sample_1", "sample_2"],
            ["sample_3", "sample_6"],
            ["sample_4", "sample_5"],
        ]

        observed = processing(distances, samples, MAX_CLUSTER_SIZE=2)
        self.assertListListEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
