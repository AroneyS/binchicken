#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from binchicken.workflow.scripts.evaluate import evaluate

SINGLEM_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"]
TARGET_COLUMNS=SINGLEM_COLUMNS+["target"]
APPRAISE_COLUMNS=SINGLEM_COLUMNS+["found_in"]
CLUSTER_COLUMNS=["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"]
EDGE_COLUMNS=["style", "cluster_size", "samples", "target_ids"]

OUTPUT_COLUMNS={
    "coassembly": str,
    "gene": str,
    "sequence": str,
    "genome": str,
    "target": str,
    "found_in": str,
    "source_samples": str,
    "source_num_hits": int,
    "source_coverage": float,
    "taxonomy": str,
    }
SUMMARY_COLUMNS = {
    "coassembly": str,
    "statistic": str,
    "within": str,
    "match": int,
    "nonmatch": int,
    "total": int,
    "match_percent": float,
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False, check_row_order=False)

    def test_evaluate_script(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 1, 2.0, "Root; old", 11],
            ["S3.1", "sample_2.1", "CCC", 1, 2.0, "Root; old", 11],
            ["S3.1", "sample_1.1", "DDD", 1, 2.0, "Root; old", 12],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_1", 5, 10.0, "Root"],
            ["coassembly_0", "S3.1", "CCC", "genome_1_transcripts", None, None, "sample_1,sample_2", 2, 4.0, "Root"],
            ["coassembly_0", "S3.1", "DDD", "genome_1_transcripts", None, None, "sample_1", 1, 2.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 1, 4, 5, 20.00],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 2, 3, 5, 40.00],
            ["coassembly_0", "novel_sequences", "recovery", 1, 4, 5, 20.00],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_all_targets(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "novel_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_none_targets(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", None, "10", None, "sample_1,sample_2", 10, 20.0, "Root; old"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 0, 1, 1, 0],
            ["coassembly_0", "taxonomy", "targets", 0, 1, 1, 0],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "novel_sequences", "recovery", 1, 0, 1, 100],
            ["coassembly_0", "bins", "recovery", 0, 2, 2, 0],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_empty_recovery(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 1, 2.0, "Root; old", 11],
            ["S3.1", "sample_2.1", "CCC", 1, 2.0, "Root; old", 11],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            }

        expected_matches = pl.DataFrame([
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_nontarget_wrong_sample(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_3.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_0", "novel_sequences", "recovery", 2, 1, 3, 66.67],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_samples_with_same_sequence(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_2.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_1,sample_2", 10, 20.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_multiple_coassemblies(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_3.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_1.1", "EEE", 1, 2.0, "Root; old", 12],
            ["S3.1", "sample_3.1", "EEE", 1, 2.0, "Root; old", 12],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_3.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_3", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
            ["match", 2, "sample_1.1,sample_3.1", "11"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "EEE", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_1", 5, 10.0, "Root"],
            ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "sample_1,sample_3", 10, 20.0, "Root"],
            ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "sample_3", 5, 10.0, "Root"],
            ["coassembly_0", "S3.1", "EEE", None, None, None, "sample_1", 1, 2.0, "Root; old"],
            ["coassembly_1", "S3.1", "EEE", "genome_1_transcripts", None, None, "sample_1,sample_3", 2, 4.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
            ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
            ["coassembly_1", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_1", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_1", "nontarget_bin_sequences", "recovery", 1, 3, 4, 25.0],
            ["coassembly_1", "nontarget_unbin_sequences", "recovery", 1, 3, 4, 25.0],
            ["coassembly_1", "novel_sequences", "recovery", 1, 3, 4, 25.0],
            ["coassembly_1", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_subset_coassembly(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_3.1", "CCC", 5, 10.0, "Root; old", 11],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_3.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2,sample_3", 3, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
            ["match", 2, "sample_1.1,sample_3.1", "11"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "CCD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_1", 5, 10.0, "Root"],
            ["coassembly_0", "S3.1", "CCC", "genome_1_transcripts", "11", None, "sample_1,sample_3", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "sample_3", 5, 10.0, "Root"],
            ["coassembly_1", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_1", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_1", 5, 10.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
            ["coassembly_0", "S3.1", "CCD", "genome_1_transcripts", None, None, None, None, None, "Root"],
            ["coassembly_1", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 2, 0, 2, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 2, 4, 6, 33.33],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 6, 6, 0],
            ["coassembly_0", "novel_sequences", "recovery", 2, 4, 6, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
            ["coassembly_1", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_1", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_1", "nontarget_bin_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "nontarget_unbin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_duplicate_sequences(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_3.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_2.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_1.1", "YYY", 5, 10.0, "Root; old", 98],
            ["S3.1", "sample_2.1", "YYY", 5, 10.0, "Root; old", 98],
            ["S3.1", "sample_2.1", "ZZZ", 5, 10.0, "Root; old", 99],
            ["S3.1", "sample_3.1", "ZZZ", 5, 10.0, "Root; old", 99],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_3.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_1.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2,sample_3", 3, 2, 3000, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_2", 2, 2, 2000, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_3.1", "10"],
            ["match", 2, "sample_1.1,sample_2.1", "11,98"],
            ["match", 2, "sample_2.1,sample_3.1", "99"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_0_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            "coassembly_1-genome_2": "genome_2.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_0_transcripts", "10", None, "sample_1,sample_3", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_3", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "sample_3", 5, 10.0, "Root"],
            ["coassembly_0", "S3.1", "CCC", None, "11", None, "sample_1,sample_2", 10, 20.0, "Root; old"],
            ["coassembly_0", "S3.1", "YYY", None, "98", None, "sample_1,sample_2", 10, 20.0, "Root; old"],
            ["coassembly_0", "S3.1", "ZZZ", None, "99", None, "sample_2,sample_3", 10, 20.0, "Root; old"],
            ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "sample_1", 5, 10.0, "Root"],
            ["coassembly_1", "S3.1", "YYY", None, "98", None, "sample_1,sample_2", 10, 20.0, "Root; old"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, None, None, None, "Root"],
            ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, None, None, None, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 3, 4, 25],
            ["coassembly_0", "taxonomy", "targets", 1, 1, 2, 50],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 1, 3, 4, 25],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 0, 4, 4, 0],
            ["coassembly_0", "novel_sequences", "recovery", 1, 3, 4, 25],
            ["coassembly_0", "bins", "recovery", 2, 0, 2, 100],
            ["coassembly_1", "sequences", "targets", 1, 1, 2, 50],
            ["coassembly_1", "taxonomy", "targets", 1, 1, 2, 50],
            ["coassembly_1", "nontarget_bin_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "nontarget_unbin_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "bins", "recovery", 1, 2, 3, 33.33],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_other_sample_target(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 1, 2.0, "Root; old", 11],
            ["S3.1", "sample_3.1", "CCC", 1, 2.0, "Root; old", 11],
        ], orient="row", schema=TARGET_COLUMNS)
        binned = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], orient="row", schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["match", 2, "sample_1.1,sample_2.1", "10"],
        ], orient="row", schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
        ], orient="row", schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "sample_1,sample_2", 10, 20.0, "Root"],
            ["coassembly_0", "S3.1", "CCC", "genome_1_transcripts", None, None, "sample_1", 1, 2.0, "Root"],
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
        ], orient="row", schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_bin_sequences", "recovery", 0, 2, 2, 0],
            ["coassembly_0", "nontarget_unbin_sequences", "recovery", 1, 1, 2, 50.0],
            ["coassembly_0", "novel_sequences", "recovery", 0, 2, 2, 0],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], orient="row", schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, binned, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)


if __name__ == '__main__':
    unittest.main()
