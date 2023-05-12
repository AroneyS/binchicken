#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.no_genomes import processing

READ_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": int,
    "taxonomy": str,
    }

APPRAISE_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": int,
    "taxonomy": str,
    "found_in": str,
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def test_no_genomes(self):
        reads = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=READ_COLUMNS)

        expected_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(reads)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_no_genomes_empty_input(self):
        reads = pl.DataFrame([
        ], schema=READ_COLUMNS)

        expected_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(reads)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_no_genomes_remove_EIF(self):
        reads = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
            ["S3.18.EIF_2_alpha", "sample_1", "NNN", 5, 10, "Root"],
        ], schema=READ_COLUMNS)

        expected_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(reads)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)


if __name__ == '__main__':
    unittest.main()
