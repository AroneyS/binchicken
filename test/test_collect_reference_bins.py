#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from unittest.mock import patch
from bird_tool_utils import in_tempdir
from binchicken.workflow.scripts.collect_reference_bins import pipeline, main

APPRAISE_COLUMNS=["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]

class Tests(unittest.TestCase):
    def test_collect_reference_bins(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set(["genome_1", "genome_2"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_underscore(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1_1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1_1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set(["genome_1", "genome_2"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_suffix_R(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1_R1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1_R1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set(["genome_1", "genome_2"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1_R1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_no_hits(self):
        appraise_binned = pl.DataFrame([
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set()
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_no_sequences(self):
        appraise_binned = pl.DataFrame([
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set()
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_low_hits(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 1, 2, "Root", "genome_1_protein,genome_2_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 10, 20, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set()
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_trimmed(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein,genome_2_protein"],
            ["S3.2", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.3", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.4", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.5", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.6", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.7", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.8", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.9", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.10", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set(["genome_1"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_multiple_samples(self):
        appraise_binned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", "genome_1_protein"],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", "genome_2_protein"],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", "genome_2_protein"],
        ], orient="row", schema=APPRAISE_COLUMNS)
        appraise_unbinned = pl.DataFrame([
            ["S3.1", "sample_1", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_2", "AAA", 5, 10, "Root", ""],
            ["S3.1", "sample_3", "AAA", 5, 10, "Root", ""],
        ], orient="row", schema=APPRAISE_COLUMNS)

        expected = set(["genome_1"])
        observed = pipeline(appraise_binned, appraise_unbinned, "sample_1")
        self.assertEqual(expected, observed)

    def test_collect_reference_bins_main_from_files(self):
        # Tables are scanned lazily and filtered to the sample; numeric sample names must still match
        with in_tempdir():
            pl.DataFrame([
                ["S3.1", "1234", "AAA", 5, 10, "Root", "genome_1_protein"],
                ["S3.1", "5678", "AAA", 5, 10, "Root", "genome_2_protein"],
            ], orient="row", schema=APPRAISE_COLUMNS).write_csv("binned.tsv", separator="\t")
            pl.DataFrame([
                ["S3.1", "1234", "CCC", 5, 10, "Root", ""],
                ["S3.1", "5678", "CCC", 5, 10, "Root", ""],
            ], orient="row", schema=APPRAISE_COLUMNS).write_csv("unbinned.tsv", separator="\t")
            for name in ["genome_1", "genome_2"]:
                with open(f"{name}.fna", "w") as f:
                    f.write(f">contig\nACGT\n")
            with open("genomes.tsv", "w") as f:
                f.write("genome_1\tgenome_1.fna\ngenome_2\tgenome_2.fna\n")

            argv = ["collect_reference_bins.py", "--appraise-binned", "binned.tsv", "--appraise-unbinned", "unbinned.tsv",
                "--genomes", "genomes.tsv", "--sample", "1234", "--read", "1234", "--output", "reference.fna"]
            with patch("sys.argv", argv):
                main()

            with open("reference.fna") as f:
                self.assertEqual(">genome_1~contig\nACGT\n", f.read())


if __name__ == '__main__':
    unittest.main()
