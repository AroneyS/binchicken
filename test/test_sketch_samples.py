#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from binchicken.workflow.scripts.sketch_samples import processing
from sourmash import load_file_as_signatures
from bird_tool_utils import in_tempdir

OTU_TABLE_COLUMNS = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

class Tests(unittest.TestCase):
    def test_sketch_samples(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"], # 3
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5

                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", threads=1)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)

    def test_sketch_samples_taxa_of_interest(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 2

                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 2

                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 3
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 4
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 5
                ["S3.1", "sample_3", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Archaea;p__Nanoarchaeota"], # 2

                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 4
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 5
                ["S3.1", "sample_4", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Archaea;p__Nanoarchaeota"], # 2
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", taxa_of_interest="p__Actinobacteriota", threads=1)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)

    def test_sketch_samples_exclude_taxa(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 2

                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 2

                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 3
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 4
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 5
                ["S3.1", "sample_3", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Archaea;p__Nanoarchaeota"], # 2

                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 1
                ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 4
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root;d__Bacteria;p__Actinobacteriota"], # 5
                ["S3.1", "sample_4", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root;d__Archaea;p__Nanoarchaeota"], # 2
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", exclude_taxa="d__Archaea", threads=1)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)

    def test_sketch_samples_exact_groups(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"], # 3
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5

                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", threads=1, samples_per_group=4)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)

    def test_sketch_samples_small_groups(self):
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_1", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_2", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"], # 2

                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"], # 3
                ["S3.1", "sample_3", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_3", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5

                ["S3.1", "sample_4", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"], # 1
                ["S3.1", "sample_4", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 4
                ["S3.1", "sample_4", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 5, 10, "Root"], # 5
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            expected_names = [
                "sample_1",
                "sample_2",
                "sample_3",
                "sample_4",
            ]

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", threads=1, samples_per_group=2)
            signatures = [s for s in load_file_as_signatures(signatures_path)]
            observed_names = [s.name for s in signatures]
            self.assertEqual(sorted(expected_names), sorted(observed_names))

            sample_1_sig = signatures[observed_names.index("sample_1")]
            sample_2_sig = signatures[observed_names.index("sample_2")]
            sample_3_sig = signatures[observed_names.index("sample_3")]
            sample_4_sig = signatures[observed_names.index("sample_4")]

            self.assertEqual(sample_1_sig.jaccard(sample_2_sig), 1.0)
            self.assertEqual(sample_1_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_1_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_2_sig.jaccard(sample_3_sig), 0.0)
            self.assertEqual(sample_2_sig.jaccard(sample_4_sig), 0.25)
            self.assertEqual(sample_3_sig.jaccard(sample_4_sig), 0.5)

    def test_sketch_samples_found_in_header(self):
        # Appraise unbinned tables carry a found_in column; the header row must not become a sample
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root", ""],
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root", ""],
            ], orient="row", schema={**OTU_TABLE_COLUMNS, "found_in": str})
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            signatures_path = processing(unbinned_path, output_path="./signatures.sig", threads=1)
            observed_names = [s.name for s in load_file_as_signatures(signatures_path)]
            self.assertEqual(["sample_1", "sample_2"], sorted(observed_names))

    def test_sketch_samples_taxa_last_row_filtered(self):
        # The final row is filtered out, which previously dropped the last kept sample
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root; d__Bacteria"],
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root; d__Bacteria"],
                ["S3.1", "sample_3", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root; d__Archaea"],
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned_path = "unbinned.otu_table.tsv"
            unbinned.write_csv(unbinned_path, separator="\t")

            for kwargs in [{"taxa_of_interest": "d__Bacteria"}, {"exclude_taxa": "d__Archaea"}]:
                signatures_path = processing(unbinned_path, output_path="./signatures.sig", threads=1, **kwargs)
                observed_names = [s.name for s in load_file_as_signatures(signatures_path)]
                self.assertEqual(["sample_1", "sample_2"], sorted(observed_names))

    def test_sketch_samples_stale_thread_files(self):
        # Leftover per-thread signature files from a previous run must not be merged into the output
        with in_tempdir():
            stale_unbinned = pl.DataFrame([
                ["S3.1", f"stale_{i}", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"]
                for i in range(3)
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            stale_unbinned.write_csv("stale.otu_table.tsv", separator="\t")
            processing("stale.otu_table.tsv", output_path="./signatures.sig", threads=1, samples_per_group=1)

            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"],
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned.write_csv("unbinned.otu_table.tsv", separator="\t")

            signatures_path = processing("unbinned.otu_table.tsv", output_path="./signatures.sig", threads=1, samples_per_group=1)
            observed_names = [s.name for s in load_file_as_signatures(signatures_path)]
            self.assertEqual(["sample_1"], observed_names)

    def test_sketch_samples_non_contiguous_sample(self):
        # A sample whose rows are split by another sample would be sketched twice from partial sequences
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root"],
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root"],
                ["S3.1", "sample_1", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root"],
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned.write_csv("unbinned.otu_table.tsv", separator="\t")

            with self.assertRaisesRegex(ValueError, "sample_1"):
                processing("unbinned.otu_table.tsv", output_path="./signatures.sig", threads=1)

    def test_sketch_samples_non_contiguous_filtered_between(self):
        # Rows of another sample that are all filtered out do not break contiguity
        with in_tempdir():
            unbinned = pl.DataFrame([
                ["S3.1", "sample_1", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 5, 10, "Root; d__Bacteria"],
                ["S3.1", "sample_2", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 5, 10, "Root; d__Archaea"],
                ["S3.1", "sample_1", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 5, 10, "Root; d__Bacteria"],
            ], orient="row", schema=OTU_TABLE_COLUMNS)
            unbinned.write_csv("unbinned.otu_table.tsv", separator="\t")

            signatures_path = processing("unbinned.otu_table.tsv", output_path="./signatures.sig", threads=1, exclude_taxa="d__Archaea")
            observed_names = [s.name for s in load_file_as_signatures(signatures_path)]
            self.assertEqual(["sample_1"], observed_names)


if __name__ == '__main__':
    unittest.main()
