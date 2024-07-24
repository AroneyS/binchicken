#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from binchicken.workflow.scripts.abundance_weighting import pipeline

APPRAISE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
}

WEIGHTING_COLUMNS = {
    "gene": str,
    "sequence": str,
    "weight": float,
}

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False, check_row_order=False)

    def test_abundance_weighting(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_other_genes(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", ""],
            ["S3.2", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.2", "sample_2", "AAA", 9, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
            ["S3.2", "AAA", 0.5],
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_no_binned(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        expected = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
            ["S3.1", "AAB", 0.5],
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_no_unbinned(self):
        unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = pl.DataFrame([
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_empty_input(self):
        unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        expected = pl.DataFrame([
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_specific_samples(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_3", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        samples = ["sample_1", "sample_2"]

        expected = pl.DataFrame([
            ["S3.1", "AAA", 0.5],
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned, samples)
        self.assertDataFrameEqual(expected, observed)

    def test_abundance_weighting_specific_samples_more(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_3", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_4", "AAA", 9, 10, "Root", ""],
            ["S3.1", "sample_4", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_4", "AAC", 5, 20, "Root", ""],
        ], schema=APPRAISE_COLUMNS)
        samples = ["sample_3", "sample_4"]

        expected = pl.DataFrame([
            ["S3.1", "AAA", 0.375],
            ["S3.1", "AAB", 0.375],
        ], schema=WEIGHTING_COLUMNS)

        observed = pipeline(unbinned, binned, samples)
        self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
