#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.target_elusive import pipeline

APPRAISE_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

TARGETS_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "target"]
EDGES_COLUMNS={
    "taxa_group": str,
    "weight": int,
    "target_ids": str,
    "sample1": str,
    "sample2": str
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False, check_row_order=False)

    def test_target_elusive(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["Root", 1, "0", "sample_1", "sample_2"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_empty_input(self):
        unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_low_coverage(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", "0"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_two_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_multiple_genes(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_three_sample_mix(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", "2"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
            ["Root", 1, "0", "sample_1", "sample_3"],
            ["Root", 1, "0", "sample_2", "sample_3"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_target_taxa(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["p__Planctomycetota", 1, "0", "sample_1", "sample_2"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, TAXA_OF_INTEREST="p__Planctomycetota")
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)


if __name__ == '__main__':
    unittest.main()
