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
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", ""],
        ], columns=APPRAISE_COLUMNS)

        expected_targets = pd.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10, "Root", "0"],
            ["S3.1", "sample_2.1", "AAA", 5, 10, "Root", "0"],
        ], columns=TARGETS_COLUMNS).astype({"target": object})
        expected_edges = pd.DataFrame([
            ["Root", 1, "0", "sample_1.1", "sample_2.1"],
        ], columns=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertEdgesDfEqual(expected_edges, observed_edges)


if __name__ == '__main__':
    unittest.main()
