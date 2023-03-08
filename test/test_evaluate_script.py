#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.evaluate import evaluate

SINGLEM_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"]
TARGET_COLUMNS=SINGLEM_COLUMNS+["target"]
APPRAISE_COLUMNS=SINGLEM_COLUMNS+["found_in"]
CLUSTER_COLUMNS=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples", "coassembly"]
EDGE_COLUMNS=["taxa_group", "weight", "target_ids", "sample1", "sample2"]

OUTPUT_COLUMNS={
    "coassembly": str,
    "gene": str,
    "sequence": str,
    "genome": str,
    "target": str,
    "found_in": str,
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
        assert_frame_equal(a, b, check_dtype=False, check_row_order=False)

    def test_evaluate_script(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_all_targets(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "novel_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_none_targets(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", None, "10", None, "Root; old"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 0, 1, 1, 0],
            ["coassembly_0", "taxonomy", "targets", 0, 1, 1, 0],
            ["coassembly_0", "nontarget_sequences", "recovery", 0, 1, 1, 0],
            ["coassembly_0", "novel_sequences", "recovery", 1, 0, 1, 100],
            ["coassembly_0", "bins", "recovery", 0, 2, 2, 0],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_empty_recovery(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            }

        expected_matches = pl.DataFrame([
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_nontarget_wrong_sample(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_3.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 0, 3, 3, 0],
            ["coassembly_0", "novel_sequences", "recovery", 2, 1, 3, 66.67],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_samples_with_same_sequence(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_2.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_multiple_coassemblies(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_3.1", "CCC", 5, 10.0, "Root; old", 11],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_3.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_3", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
            ["Root", 1, "11", "sample_1.1", "sample_3.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
            ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
            ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
            ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
            ["coassembly_1", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_1", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_1", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)

    def test_evaluate_script_subset_coassembly(self):
        targets = pl.DataFrame([
            ["S3.1", "sample_1.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_2.1", "AAA", 5, 10.0, "Root; old", 10],
            ["S3.1", "sample_1.1", "CCC", 5, 10.0, "Root; old", 11],
            ["S3.1", "sample_3.1", "CCC", 5, 10.0, "Root; old", 11],
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_1.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_3.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2,sample_3", 3, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_2", 2, 1, 1, 100, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_2.1"],
            ["Root", 1, "11", "sample_1.1", "sample_3.1"],
        ], schema=EDGE_COLUMNS)
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
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
            ["coassembly_0", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
            ["coassembly_0", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
            ["coassembly_1", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_1", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
            ["coassembly_0", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
            ["coassembly_1", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 2, 0, 2, 100],
            ["coassembly_0", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_0", "nontarget_sequences", "recovery", 2, 4, 6, 33.33],
            ["coassembly_0", "novel_sequences", "recovery", 2, 4, 6, 33.33],
            ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
            ["coassembly_1", "sequences", "targets", 1, 0, 1, 100],
            ["coassembly_1", "taxonomy", "targets", 1, 0, 1, 100],
            ["coassembly_1", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "bins", "recovery", 1, 1, 2, 50],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
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
        ], schema=TARGET_COLUMNS)
        nontargets = pl.DataFrame([
            ["S3.1", "sample_3.1", "BBB", 5, 10.0, "Root; old", "oldgenome_1"],
            ["S3.1", "sample_1.1", "DDD", 5, 10.0, "Root; old", "oldgenome_2"],
        ], schema=APPRAISE_COLUMNS)
        clusters = pl.DataFrame([
            ["sample_1,sample_2,sample_3", 3, 2, 2, 3000, "sample_1,sample_2,sample_3", "coassembly_0"],
            ["sample_1,sample_2", 2, 2, 2, 2000, "sample_1,sample_2,sample_3", "coassembly_1"],
        ], schema=CLUSTER_COLUMNS)
        edges = pl.DataFrame([
            ["Root", 1, "10", "sample_1.1", "sample_3.1"],
            ["Root", 2, "11,98", "sample_1.1", "sample_2.1"],
            ["Root", 1, "99", "sample_2.1", "sample_3.1"],
        ], schema=EDGE_COLUMNS)
        recovered = pl.DataFrame([
            ["S3.1", "coassembly_0-genome_0_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAA", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "AAB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_0-genome_1_transcripts", "BBB", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCC", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "CCD", 1, 2.0, "Root"],
            ["S3.1", "coassembly_1-genome_1_transcripts", "DDD", 1, 2.0, "Root"],
        ], schema=SINGLEM_COLUMNS)
        recovered_bins = {
            "coassembly_0-genome_0": "genome_0.fna",
            "coassembly_0-genome_1": "genome_1.fna",
            "coassembly_1-genome_0": "genome_0.fna",
            "coassembly_1-genome_1": "genome_1.fna",
            "coassembly_1-genome_2": "genome_2.fna",
            }

        expected_matches = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAA", "genome_0_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
            ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
            ["coassembly_0", "S3.1", "CCC", None, "11", None, "Root; old"],
            ["coassembly_0", "S3.1", "YYY", None, "98", None, "Root; old"],
            ["coassembly_0", "S3.1", "ZZZ", None, "99", None, "Root; old"],
            ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
            ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
            ["coassembly_1", "S3.1", "YYY", None, "98", None, "Root; old"],
        ], schema=OUTPUT_COLUMNS)
        expected_unmatched = pl.DataFrame([
            ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
            ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
        ], schema=OUTPUT_COLUMNS)
        expected_summary = pl.DataFrame([
            ["coassembly_0", "sequences", "targets", 1, 3, 4, 25],
            ["coassembly_0", "taxonomy", "targets", 1, 1, 2, 50],
            ["coassembly_0", "nontarget_sequences", "recovery", 1, 3, 4, 25],
            ["coassembly_0", "novel_sequences", "recovery", 1, 3, 4, 25],
            ["coassembly_0", "bins", "recovery", 2, 0, 2, 100],
            ["coassembly_1", "sequences", "targets", 1, 1, 2, 50],
            ["coassembly_1", "taxonomy", "targets", 1, 1, 2, 50],
            ["coassembly_1", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
            ["coassembly_1", "bins", "recovery", 1, 2, 3, 33.33],
        ], schema=SUMMARY_COLUMNS)

        observed_matches, observed_unmatched, observed_summary = evaluate(targets, nontargets, clusters, edges, recovered, recovered_bins)
        self.assertDataFrameEqual(expected_matches, observed_matches)
        self.assertDataFrameEqual(expected_unmatched, observed_unmatched)
        self.assertDataFrameEqual(expected_summary, observed_summary)


if __name__ == '__main__':
    unittest.main()
