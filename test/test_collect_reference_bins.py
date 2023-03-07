#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from ibis.workflow.scripts.collect_reference_bins import pipeline

APPRAISE_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

class Tests(unittest.TestCase):
    def test_collect_reference_bins(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
        ], schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = set(["genome_1", "genome_2"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_no_hits(self):
        appraise_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = set()
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_low_hits(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 1, 2, "Root", "genome_1_protein,genome_2_protein"],
        ], schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 10, 20, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = set()
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_trimmed(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
            ["S3.2", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.3", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.4", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.5", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.6", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.7", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.8", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.9", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.10", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
        ], schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = set(["genome_1"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_multiple_samples(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", "genome_2_protein"],
            ["S3.1", "sample_3.1", "AAA", 5, 10, "Root", "genome_2_protein"],
        ], schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3.1", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected = set(["genome_1"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
