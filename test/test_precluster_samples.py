#!/usr/bin/env python3

import unittest
import numpy as np
from ibis.workflow.scripts.precluster_samples import processing

class Tests(unittest.TestCase):
    def assertListListEqual(self, a, b):
        a_sets = [set(x) for x in a]
        b_sets = [set(x) for x in b]
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
            ["sample_3", "sample_4"]
        ]

        observed = processing(distances, samples, MAX_CLUSTER_SIZE=2)
        self.assertListListEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
