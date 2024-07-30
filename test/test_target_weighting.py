#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from binchicken.workflow.scripts.target_weighting import pipeline

TARGET_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
    "target": str,
}

WEIGHTING_COLUMNS = {
    "gene": str,
    "sequence": str,
    "weight": float,
}

TARGET_WEIGHTING_COLUMNS = {
    "target": str,
    "weight": float,
}

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False, check_row_order=False)

    def test_target_weighting(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", "", "0"],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", "", "0"],
            ["S3.1", "sample_2", "BBB", 9, 10, "Root", "", "1"],
        ], schema=TARGET_COLUMNS)
        weighting = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
            ["S3.2", "AAA", 0.9],
            ["S3.1", "BBB", 0.9],
        ], schema=WEIGHTING_COLUMNS)

        expected = pl.DataFrame([
            ["0", 0.5],
            ["1", 0.9],
        ], schema=TARGET_WEIGHTING_COLUMNS)

        observed = pipeline(targets, weighting)
        self.assertDataFrameEqual(expected, observed)

    def test_target_weighting_empty_input(self):
        targets = pl.DataFrame([
        ], schema=TARGET_COLUMNS)
        weighting = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
        ], schema=WEIGHTING_COLUMNS)

        expected = pl.DataFrame([
        ], schema=TARGET_WEIGHTING_COLUMNS)

        observed = pipeline(targets, weighting)
        self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
