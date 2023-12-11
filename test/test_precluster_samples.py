#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from ibis.workflow.scripts.precluster_samples import processing

OTU_TABLE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

class Tests(unittest.TestCase):
    def assertListListEqual(self, a, b):
        a_sets = [set(x) for x in a]
        b_sets = [set(x) for x in b]
        self.assertEqual(sorted(a_sets), sorted(b_sets))

    def test_precluster_samples(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
            ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

            ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
            ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

            ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"], # 3
            ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
            ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5

            ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
            ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
            ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5
        ], schema=OTU_TABLE_COLUMNS)

        expected = [
            ["sample_1", "sample_2"],
            ["sample_3", "sample_4"]
        ]

        observed = processing(unbinned, MAX_CLUSTER_SIZE=2)
        self.assertListListEqual(expected, observed)

    def test_precluster_samples(self):
        sequences = [
            "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA",
            "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC",
            "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA",
            "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG",
            "ATCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA",
            "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
        ]
        # large group of 99 samples with 3 in common
        large = [
            ["S3.1", f"sample_{i}", s, 5, 10, "Root"] for s in sequences[0:3] for i in range(1, 100)
        ]
        # another group of 99 samples with 3 in common (2 overlap above)
        overlap = [
            ["S3.1", f"sample_{i}", s, 5, 10, "Root"] for s in sequences[1:4] for i in range(100, 199)
        ]
        # another group of 20 samples with 2 in common (independent)
        small = [
            ["S3.1", f"sample_{i}", s, 5, 10, "Root"] for s in sequences[4:6] for i in range(199, 218)
        ]

        unbinned = pl.DataFrame(large + overlap + small, schema=OTU_TABLE_COLUMNS)

        observed = processing(unbinned, MAX_CLUSTER_SIZE=100)
        largest_observed = max(observed, key=len)
        smallest_observed = min(observed, key=len)
        self.assertTrue(len(largest_observed) <= 1000)
        self.assertTrue(len(smallest_observed) >= 10)


if __name__ == '__main__':
    unittest.main()
