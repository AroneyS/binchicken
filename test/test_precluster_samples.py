#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.precluster_samples import processing
from bird_tool_utils import in_tempdir

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
        with in_tempdir():
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

            observed = processing(
                unbinned,
                os.path.join(os.getcwd(), "test_precluster_samples"),
                MAX_CLUSTER_SIZE=2,
                )
            self.assertListListEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
