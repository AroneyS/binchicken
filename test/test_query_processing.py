#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.query_processing import processing

QUERY_COLUMNS = {
    "query_name": str,
    "query_sequence": str,
    "divergence": int,
    "num_hits": int,
    "coverage": int,
    "sample": str,
    "marker": str,
    "hit_sequence": str,
    "taxonomy": str,
    }
PIPE_COLUMNS = {
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

    def test_query_processing(self):
        query = pl.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_empty_input(self):
        query = pl.DataFrame([
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_no_binned(self):
        query = pl.DataFrame([
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_no_unbinned(self):
        query = pl.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_extra_pipe_entries(self):
        query = pl.DataFrame([
            ["sample_1", "AAA", 1, 20, 50, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_extra_entries(self):
        query = pl.DataFrame([
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAB", "Root"],
            ["sample_1", "AAA", 10, 5, 10, "genome_2", "S3.1", "AAA", "Root; Eukaryote"],
            ["sample_1", "BAB", 1, 5, 10, "genome_1", "S3.1", "BAB", "Root"],
            ["sample_1", "BAA", 1, 5, 10, "genome_2", "S3.1", "BAA", "Root; Eukaryote"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
            ["S3.1", "sample_1", "BAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "BAB", 5, 10, "Root", "genome_1"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_multiple_genomes(self):
        query = pl.DataFrame([
            ["sample_1", "AAA", 1, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAA", 1, 5, 10, "genome_2", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 10, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1,genome_2"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)

    def test_query_processing_sequence_identity(self):
        query = pl.DataFrame([
            ["sample_1", "AAA", 3, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
            ["sample_1", "AAB", 4, 5, 10, "genome_1", "S3.1", "AAA", "Root"],
        ], schema=QUERY_COLUMNS)
        pipe = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root"],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root"],
        ], schema=PIPE_COLUMNS)

        expected_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1"],
        ], schema=APPRAISE_COLUMNS)
        expected_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", None],
        ], schema=APPRAISE_COLUMNS)

        observed_binned, observed_unbinned = processing(query, pipe, SEQUENCE_IDENTITY = 0.94)
        self.assertDataFrameEqual(expected_binned, observed_binned)
        self.assertDataFrameEqual(expected_unbinned, observed_unbinned)


if __name__ == '__main__':
    unittest.main()
