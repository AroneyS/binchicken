#!/usr/bin/env python3

import unittest
import pandas as pd
from cockatoo.workflow.scripts.target_elusive import pipeline

APPRAISE_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

EDGES_COLUMNS=["taxa_group", "weight", "target_ids", "sample1", "sample2"]
TARGETS_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "target"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def combine_samples(self, df):
        df["samples"] = df[["sample1", "sample2"]].agg(sorted, axis=1)
        df.drop(columns=["sample1", "sample2"], inplace=True)
        return df

    def assertEdgesDfEqual(self, a, b):
        pd.testing.assert_frame_equal(self.combine_samples(a), self.combine_samples(b))

    def test_target_elusive(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["Root", 1, "0", "sample_1", "sample_2"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_empty_input(self):
        unbinned = pd.DataFrame([
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
        ], columns=TARGETS_COLUMNS)
        expected_edges = pd.DataFrame([
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_low_coverage(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 5, "Root", "0"],
            ["S3.1", "sample_2", "AAA", 5, 4, "Root", "0"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
        ], columns=EDGES_COLUMNS).astype({"weight": int})

        observed_targets, observed_edges = pipeline(unbinned)
        observed_edges.index = []
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_two_targets(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_multiple_genes(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_three_sample_mix(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_3", "AAC", 5, 10, "Root", "2"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["Root", 2, "0,1", "sample_1", "sample_2"],
            ["Root", 1, "0", "sample_1", "sample_3"],
            ["Root", 1, "0", "sample_2", "sample_3"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)

    def test_target_elusive_target_taxa(self):
        unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", "1"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["p__Planctomycetota", 1, "0", "sample_1", "sample_2"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, TAXA_OF_INTEREST="p__Planctomycetota")
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)


if __name__ == '__main__':
    unittest.main()
