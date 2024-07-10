#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from binchicken.workflow.scripts.sketch_samples import processing
from sourmash import load_file_as_signatures
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
    def test_sketch_samples(self):
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

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned, output_path="./signatures.sig", threads=1)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)


if __name__ == '__main__':
    unittest.main()
