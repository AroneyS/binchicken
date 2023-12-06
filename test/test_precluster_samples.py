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

CLUSTER_COLUMNS = {
    "cluster": pl.List(str),
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def test_precluster_samples(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"],
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"],
                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"],
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"],
                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"],
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"],
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"],
                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"],
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"],
            ], schema=OTU_TABLE_COLUMNS)

            expected = pl.DataFrame([
                [["sample_1", "sample_2"]]
            ], schema=CLUSTER_COLUMNS)

            observed = processing(unbinned, os.path.join(os.getcwd(), "test_precluster_samples"))
            self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
