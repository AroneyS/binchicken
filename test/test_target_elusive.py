#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "8"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.target_elusive import pipeline

APPRAISE_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "found_in": str,
}

TARGETS_COLUMNS={
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    "target": str,
}
EDGES_COLUMNS={
    "style": str,
    "cluster_size": int,
    "samples": str,
    "target_ids": str,
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
            ["match", 2, "sample_1,sample_2", "0"],
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
            ["match", 2, "sample_1,sample_2", "0,1"],
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
            ["match", 2, "sample_1,sample_2", "0,1"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_multiple_genes_same_sequence(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.2", "sample_2", "AAA", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_1", "AAA", 5, 10, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "0"],
            ["S3.2", "sample_2", "AAA", 5, 10, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0,1"],
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
            ["match", 2, "sample_1,sample_2", "0,1"],
            ["match", 2, "sample_1,sample_3", "0"],
            ["match", 2, "sample_2,sample_3", "0"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_three_way_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 3, 6, "Root", ""],
            ["S3.1", "sample_1", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAA", 3, 6, "Root", ""],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 3, 6, "Root", "0"],
            ["S3.1", "sample_1", "AAB", 2, 4, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 3, 6, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
            ["pool", 3, "sample_1,sample_2,sample_3", "0,1"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=3)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_four_way_targets(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_1", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_2", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_2", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_3", "AAC", 1, 3, "Root", ""],
            ["S3.1", "sample_4", "AAB", 2, 4, "Root", ""],
            ["S3.1", "sample_4", "AAC", 1, 3, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_1", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_2", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_2", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_2", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_3", "AAA", 2, 4, "Root", "0"],
            ["S3.1", "sample_3", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_3", "AAC", 1, 3, "Root", "1"],
            ["S3.1", "sample_4", "AAB", 2, 4, "Root", "2"],
            ["S3.1", "sample_4", "AAC", 1, 3, "Root", "1"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["pool", 3, "sample_1,sample_2,sample_3", "0"],
            ["pool", 3, "sample_2,sample_3,sample_4", "2"],
            ["pool", 4, "sample_1,sample_2,sample_3,sample_4", "1"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    @unittest.skip("Benchmarking")
    def test_target_elusive_benchmarking(self):
        genes = ["S3.1", "S3.2", "S3.3", "S3.4", "S3.5", "S3.6", "S3.7", "S3.8", "S3.9", "S3.10"]
        samples = ["sample_" + str(n) for n in range(1, 51)]
        sequences = ["AAA", "AAB", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAJ"]
        unbinned = pl.DataFrame([
            [g, s, n, 2, 5.1, "Root", ""] for g in genes for s in samples for n in sequences
        ], schema=APPRAISE_COLUMNS)

        _, _ = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)
        _, _ = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=4)

        # Cross join, recursive time: 247.833s
        # Pooling time: 1.979s

    def test_target_elusive_target_taxa(self):
        unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_1", "AAB", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", ""],
            ["S3.1", "sample_2", "AAB", 5, 10, "Root", ""],
        ], schema=APPRAISE_COLUMNS)

        expected_targets = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root; d__Bacteria; p__Planctomycetota", "0"],
        ], schema=TARGETS_COLUMNS)
        expected_edges = pl.DataFrame([
            ["match", 2, "sample_1,sample_2", "0"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, TAXA_OF_INTEREST="p__Planctomycetota")
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)

    def test_target_elusive_single_assembly(self):
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
            ["match", 2, "sample_1,sample_2", "0,1"],
            ["match", 2, "sample_1,sample_3", "0"],
            ["match", 2, "sample_2,sample_3", "0"],
        ], schema=EDGES_COLUMNS)

        observed_targets, observed_edges = pipeline(unbinned, MAX_COASSEMBLY_SAMPLES=1)
        self.assertDataFrameEqual(expected_targets, observed_targets)
        self.assertDataFrameEqual(expected_edges, observed_edges)


if __name__ == '__main__':
    unittest.main()
