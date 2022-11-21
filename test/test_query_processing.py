#!/usr/bin/env python3

import unittest
import pandas as pd
from cockatoo.workflow.scripts.query_processing import processing

QUERY_COLUMNS = ["query_name", "query_sequence", "divergence", "num_hits", "coverage", "sample", "marker", "hit_sequence", "taxonomy"]
PIPE_COLUMNS = ["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"]

APPRAISE_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def test_query_processing(self):
        query = pd.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], columns=APPRAISE_COLUMNS)
        expected_unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], columns=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)


if __name__ == '__main__':
    unittest.main()
