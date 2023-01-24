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

    def test_query_processing_empty_input(self):
        query = pd.DataFrame([
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
        ], columns=APPRAISE_COLUMNS)
        expected_unbinned = pd.DataFrame([
        ], columns=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_no_binned(self):
        query = pd.DataFrame([
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
        ], columns=APPRAISE_COLUMNS).astype({"num_hits": int, "coverage": int})
        expected_unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], columns=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        observed_binned.index = []
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_no_unbinned(self):
        query = pd.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], columns=APPRAISE_COLUMNS)
        expected_unbinned = pd.DataFrame([
        ], columns=APPRAISE_COLUMNS).astype({"num_hits": int, "coverage": int})

        observed_binned, observed_unbinned = processing(query, pipe)
        observed_unbinned.index = []
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_extra_pipe_entries(self):
        query = pd.DataFrame([
            ["sample_1", "AAA", 1, 20, 50, "genome_1", "S3.1", "AAA", "Root"],
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

    def test_query_processing_extra_entries(self):
        query = pd.DataFrame([
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAB", "Root"],
            ["sample_1", "AAA", 10, 5, 10, "genome_2", "S3.1", "AAA", "Root; Eukaryote"],
            ["sample_1", "BAB", 1, 5, 10, "genome_1", "S3.1", "BAB", "Root"],
            ["sample_1", "BAA", 1, 5, 10, "genome_2", "S3.1", "BAA", "Root; Eukaryote"],
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
            ["S3.1", "sample_1", "BAB", 5, 10, "Root"],
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
            ["S3.1", "sample_1", "BAB", 5, 10, "Root", "genome_1"],
        ], columns=APPRAISE_COLUMNS)
        expected_unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], columns=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_multiple_genomes(self):
        query = pd.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAA", 1, 5, 10, "genome_2", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], columns=QUERY_COLUMNS)
        pipe = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], columns=PIPE_COLUMNS)

        expected_binned = pd.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1,genome_2"],
        ], columns=APPRAISE_COLUMNS)
        expected_unbinned = pd.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], columns=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_sequence_identity(self):
        query = pd.DataFrame([
            ["sample_1", "AAA", 3, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 4, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
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

        observed_binned, observed_unbinned = processing(query, pipe, SEQUENCE_IDENTITY = 0.94)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)


if __name__ == '__main__':
    unittest.main()
