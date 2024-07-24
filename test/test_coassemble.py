#!/usr/bin/env python3

import unittest
import os
import gzip
from bird_tool_utils import in_tempdir
import extern
import subprocess
from snakemake.io import load_configfile
import polars as pl
from polars.testing import assert_frame_equal

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

SAMPLE_READS_FORWARD = " ".join([
    os.path.join(path_to_data, "sample_1.1.fq"),
    os.path.join(path_to_data, "sample_2.1.fq"),
    os.path.join(path_to_data, "sample_3.1.fq"),
])
SAMPLE_READS_REVERSE = " ".join([
    os.path.join(path_to_data, "sample_1.2.fq"),
    os.path.join(path_to_data, "sample_2.2.fq"),
    os.path.join(path_to_data, "sample_3.2.fq"),
])
SAMPLE_READS_FORWARD_EMPTY = " ".join([SAMPLE_READS_FORWARD, os.path.join(path_to_data, "sample_4.1.fq")])
SAMPLE_READS_REVERSE_EMPTY = " ".join([SAMPLE_READS_REVERSE, os.path.join(path_to_data, "sample_4.2.fq")])
SAMPLE_READS_FORWARD_PRE = " ".join([SAMPLE_READS_FORWARD, os.path.join(path_to_data, "sample_5.1.fq")])
SAMPLE_READS_REVERSE_PRE = " ".join([SAMPLE_READS_REVERSE, os.path.join(path_to_data, "sample_5.2.fq")])

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
METAPACKAGE_EIF = os.path.join(path_to_data, "singlem_metapackage_EIF.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")

SAMPLE_READ_SIZE = os.path.join(MOCK_COASSEMBLE, "coassemble", "read_size2.csv")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_3_read.otu_table.tsv"),
    ])
SAMPLE_SINGLEM_EIF = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_EIF", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_EIF", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_EIF", "sample_3_read.otu_table.tsv"),
    ])
SAMPLE_SINGLEM_PRE = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_3_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_5_read.otu_table.tsv"),
    ])
SAMPLE_SINGLEM_WEIGHT = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_weight", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_weight", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_weight", "sample_3_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_weight", "sample_5_read.otu_table.tsv"),
    ])

GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "coassemble", "summarise", "bins_summarised.otu_table2.tsv")
GENOME_SINGLEM_EIF = os.path.join(MOCK_COASSEMBLE, "coassemble", "summarise", "bins_summarised_EIF.otu_table.tsv")

SAMPLE_QUERY_DIR = os.path.join(path_to_data, "query")
SAMPLE_QUERY = ' '.join([
    os.path.join(path_to_data, "query", "sample_1_query.otu_table.tsv"),
    os.path.join(path_to_data, "query", "sample_2_query.otu_table.tsv"),
    os.path.join(path_to_data, "query", "sample_3_query.otu_table.tsv"),
    ])
SAMPLE_QUERY_SINGLEM = ' '.join([
    os.path.join(path_to_data, "query", "sample_1_read.otu_table.tsv"),
    os.path.join(path_to_data, "query", "sample_2_read.otu_table.tsv"),
    os.path.join(path_to_data, "query", "sample_3_read.otu_table.tsv"),
    ])

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

class Tests(unittest.TestCase):
    def test_coassemble(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_EMPTY} "
                f"--reverse {SAMPLE_READS_REVERSE_EMPTY} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--assemble-unmapped "
                f"--unmapping-max-identity 99 "
                f"--unmapping-max-alignment 90 "
                f"--prodigal-meta "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            self.assertTrue("1 samples had no targets with sufficient combined coverage for coassembly prediction" in output)
            self.assertTrue("These are recorded at " in output)

            unused_samples_path = os.path.join("test", "coassemble", "target", "unused_samples.tsv")
            self.assertTrue(os.path.exists(unused_samples_path))
            expected = "\n".join(["sample_4"])
            with open(unused_samples_path) as f:
                self.assertEqual(expected, f.read())

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "4832"]),
                    ",".join(["sample_2", "3926"]),
                    ",".join(["sample_3", "3624"]),
                    ",".join(["sample_4", "604"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
                self.assertEqual(expected, f.read())

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_2",
                        "2",
                        "3",
                        "8758",
                        "sample_1,sample_2,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

            bins_reference_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertFalse(os.path.exists(bins_reference_path))

            output_bam_files = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.bam")
            self.assertFalse(os.path.exists(output_bam_files))

            coverm_working_dir = os.path.join("test", "coassemble", "mapping", "sample_1_coverm")
            self.assertFalse(os.path.exists(coverm_working_dir))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:112:HMNM5DSXX:4:9999:19126:17300" not in file)

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            test_dir = os.path.abspath("test")
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble --coassemble -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble"),
                        "-n 64 -t 64 -m 500 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "recover"),
                        "--binning-only --refinery-max-iterations 0 "
                        "-n 32 -t 32 -m 250 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

            summary_path = os.path.join("test", "coassemble", "summary.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "unmapped_size",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sample_1,sample_2",
                        "2",
                        "3",
                        "8758",
                        "8456",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_taxa_of_interest(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--taxa-of-interest \"p__Actinobacteriota\" "
                f"--aviary-speed comprehensive "
                f"--run-qc "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_3",
                        "2",
                        "2",
                        "2869",
                        "sample_1,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

            qc_sample_1F_path = os.path.join("test", "coassemble", "qc", "sample_1_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_1F_path))
            qc_sample_1R_path = os.path.join("test", "coassemble", "qc", "sample_1_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_1R_path))

            qc_sample_2F_path = os.path.join("test", "coassemble", "qc", "sample_2_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_2F_path))
            qc_sample_2R_path = os.path.join("test", "coassemble", "qc", "sample_2_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_2R_path))

            qc_sample_3F_path = os.path.join("test", "coassemble", "qc", "sample_3_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_3F_path))
            qc_sample_3R_path = os.path.join("test", "coassemble", "qc", "sample_3_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_3R_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            test_dir = os.path.abspath("test")
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble --coassemble -1",
                        os.path.join(test_dir, "coassemble", "qc", "sample_1_1.fastq.gz"),
                        os.path.join(test_dir, "coassemble", "qc", "sample_3_1.fastq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "qc", "sample_1_2.fastq.gz"),
                        os.path.join(test_dir, "coassemble", "qc", "sample_3_2.fastq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble"),
                        "-n 64 -t 64 -m 500 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(test_dir, "coassemble", "qc", "sample_1_1.fastq.gz"),
                        os.path.join(test_dir, "coassemble", "qc", "sample_3_1.fastq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "qc", "sample_1_2.fastq.gz"),
                        os.path.join(test_dir, "coassemble", "qc", "sample_3_2.fastq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "recover"),
                        "-n 32 -t 32 -m 250 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_unmap_runqc(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--assemble-unmapped "
                f"--run-qc "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            qc_sample_1F_path = os.path.join("test", "coassemble", "qc", "sample_1_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_1F_path))
            qc_sample_1R_path = os.path.join("test", "coassemble", "qc", "sample_1_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_1R_path))

            qc_sample_2F_path = os.path.join("test", "coassemble", "qc", "sample_2_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_2F_path))
            qc_sample_2R_path = os.path.join("test", "coassemble", "qc", "sample_2_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_2R_path))

            qc_sample_3F_path = os.path.join("test", "coassemble", "qc", "sample_3_1.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_3F_path))
            qc_sample_3R_path = os.path.join("test", "coassemble", "qc", "sample_3_2.fastq.gz")
            self.assertTrue(os.path.exists(qc_sample_3R_path))

            map_sample_1F_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(map_sample_1F_path))
            map_sample_1R_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.2.fq.gz")
            self.assertTrue(os.path.exists(map_sample_1R_path))

            map_sample_2F_path = os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(map_sample_2F_path))
            map_sample_2R_path = os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.2.fq.gz")
            self.assertTrue(os.path.exists(map_sample_2R_path))

            map_sample_3F_path = os.path.join("test", "coassemble", "mapping", "sample_3_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(map_sample_3F_path))
            map_sample_3R_path = os.path.join("test", "coassemble", "mapping", "sample_3_unmapped.2.fq.gz")
            self.assertTrue(os.path.exists(map_sample_3R_path))

    def test_coassemble_query_input(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-query {SAMPLE_QUERY} "
                f"--sample-singlem {SAMPLE_QUERY_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"cluster_graph\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_2",
                        "2",
                        "3",
                        "2869",
                        "sample_1,sample_2,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_query_input_taxa_of_interest(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {TWO_GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-query {SAMPLE_QUERY} "
                f"--sample-singlem {SAMPLE_QUERY_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--taxa-of-interest \"p__Actinobacteriota\" "
                f"--assemble-unmapped "
                f"--unmapping-max-identity 99 "
                f"--unmapping-max-alignment 90 "
                f"--unmapping-min-appraised 0.09 "
                f"--prodigal-meta "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "2", "3.28", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter sp013286235", "GB_GCA_013286235.1,GB_GCA_013286235.2"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter sp013286235", "GB_GCA_013286235.1"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter sp013286235", "GB_GCA_013286235.1"]),
                    ""
                ]
            )
            with open(binned_path) as f:
                self.assertEqual(expected, f.read())

            unbinned_path = os.path.join("test", "coassemble", "appraise", "unbinned.otu_table.tsv")
            self.assertTrue(os.path.exists(unbinned_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "gene",
                        "sample",
                        "sequence",
                        "num_hits",
                        "coverage",
                        "taxonomy",
                        "found_in",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_1",
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "3",
                        "4.92",
                        "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_1",
                        "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC",
                        "5",
                        "8.21",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_1",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA",
                        "1",
                        "1.64",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_1",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATC",
                        "5",
                        "8.21",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis3",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_1",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATG",
                        "5",
                        "8.21",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis4",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_2",
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "4",
                        "6.57",
                        "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_2",
                        "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC",
                        "3",
                        "4.92",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_2",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATG",
                        "5",
                        "8.21",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis4",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_3",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA",
                        "6",
                        "9.85",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2",
                        "",
                    ]),
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        "sample_3",
                        "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATG",
                        "5",
                        "8.21",
                        "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis4",
                        "",
                    ]),
                    ""
                ]
            )
            with open(unbinned_path) as f:
                self.assertEqual(expected, f.read())

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:112:HMNM5DSXX:4:9999:19126:17300" not in file)

            summary_path = os.path.join("test", "coassemble", "summary.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected1 = "\n".join(
                [
                    "\t".join(["coassembly", "samples", "length", "total_targets", "total_size", "unmapped_size",]),
                    "\t".join(["coassembly_0", "sample_1,sample_2", "2", "2", "2869", "8456",]),
                    ""
                ]
            )
            expected2 = "\n".join(
                [
                    "\t".join(["coassembly", "samples", "length", "total_targets", "total_size", "unmapped_size",]),
                    "\t".join(["coassembly_0", "sample_1,sample_3", "2", "2", "2869", "8154",]),
                    ""
                ]
            )
            with open(summary_path) as f:
                observed = f.read()
                self.assertTrue(expected1 == observed or expected2 == observed)

    def test_coassemble_exclude_coassemblies(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--exclude-coassemblies sample_1,sample_2 "
                f"--assemble-unmapped "
                f"--unmapping-max-identity 99 "
                f"--unmapping-max-alignment 90 "
                f"--prodigal-meta "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_3",
                        "2",
                        "2",
                        "8456",
                        "sample_1,sample_2,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_single_assembly(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--single-assembly "
                f"--coassembly-samples sample_1 sample_2 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"cluster_graph\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "4",
                        "4832",
                        "sample_1,sample_2,sample_3",
                        "coassembly_0"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "3",
                        "3926",
                        "sample_1,sample_2,sample_3",
                        "coassembly_1"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_no_genomes(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"cluster_graph\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "gene",
                        "sample",
                        "sequence",
                        "num_hits",
                        "coverage",
                        "taxonomy",
                        "found_in",
                    ]),
                    ""
                ]
            )
            with open(binned_path) as f:
                self.assertEqual(expected, f.read())

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_2",
                        "2",
                        "3",
                        "8758",
                        "sample_1,sample_2,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_no_mapping(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            test_dir = os.path.abspath("test")

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble --coassemble -1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble"),
                        "-n 64 -t 64 -m 500 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_3.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_3.2.fq"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "recover"),
                        "--binning-only --refinery-max-iterations 0 "
                        "-n 32 -t 32 -m 250 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_genome_trim(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {TWO_GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--assemble-unmapped "
                f"--prodigal-meta "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"--notemp finish_mapping\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            genomes_to_map_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertTrue(os.path.exists(genomes_to_map_path))
            with open(genomes_to_map_path) as f:
                lines = f.read()
                self.assertTrue("JABDGE010000038.1_13" in lines)
                self.assertTrue("TEST_CONTIG" not in lines)

    def test_coassemble_default_config(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["taxa_of_interest"], "")
            self.assertEqual(config["assemble_unmapped"], False)
            self.assertEqual(config["aviary_assemble_threads"], 64)
            self.assertEqual(config["aviary_assemble_memory"], 500)
            self.assertEqual(config["aviary_recover_threads"], 32)
            self.assertEqual(config["aviary_recover_memory"], 250)
            self.assertEqual(config["coassembly_samples"], [])
            self.assertEqual(config["assembly_strategy"], "dynamic")

    def test_coassemble_singlem_inputs(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_run_aviary(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--run-aviary "
                f"--aviary-gtdbtk-db gtdb_release "
                f"--aviary-checkm2-db CheckM2_database "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" not in output)
            self.assertTrue("aviary_assemble" not in output)
            self.assertTrue("aviary_recover" not in output)
            self.assertTrue("aviary_combine" in output)

    def test_coassemble_run_aviary_unmapped(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--assemble-unmapped "
                f"--run-aviary "
                f"--aviary-gtdbtk-db gtdb_release "
                f"--aviary-checkm2-db CheckM2_database "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" not in output)
            self.assertTrue("aviary_assemble" not in output)
            self.assertTrue("aviary_recover" not in output)
            self.assertTrue("aviary_combine" in output)

    def test_coassemble_run_aviary_unmapped_qc(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--assemble-unmapped "
                f"--run-qc "
                f"--run-aviary "
                f"--aviary-gtdbtk-db gtdb_release "
                f"--aviary-checkm2-db CheckM2_database "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" not in output)
            self.assertTrue("aviary_assemble" not in output)
            self.assertTrue("aviary_recover" not in output)
            self.assertTrue("aviary_combine" in output)

    def test_coassemble_files_of_paths(self):
        with in_tempdir():
            write_string_to_file(SAMPLE_READS_FORWARD, "sample_reads_forward")
            write_string_to_file(SAMPLE_READS_REVERSE, "sample_reads_reverse")
            write_string_to_file(GENOMES, "genomes")
            write_string_to_file(GENOME_TRANSCRIPTS, "genome_transcripts")
            write_string_to_file("sample_1\nsample_2", "coassembly_samples")
            write_string_to_file("sample_1,sample_2", "exclude_coassemblies")
            write_string_to_file("sample_2\nsample_3", "abundance_samples")

            cmd = (
                f"binchicken coassemble "
                f"--forward-list sample_reads_forward "
                f"--reverse-list sample_reads_reverse "
                f"--genomes-list genomes "
                f"--genome-transcripts-list genome_transcripts "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--coassembly-samples-list coassembly_samples "
                f"--exclude-coassemblies-list exclude_coassemblies "
                f"--abundance-weighted-samples-list abundance_samples "
                f"--assemble-unmapped "
                f"--assembly-strategy megahit "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" in output)
            self.assertTrue("singlem_summarise_genomes" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" in output)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["genomes"], {os.path.splitext(os.path.basename(g))[0]: g for g in GENOMES.split(" ")})
            self.assertEqual(config["coassembly_samples"], ["sample_1", "sample_2"])
            self.assertEqual(config["exclude_coassemblies"], ["sample_1,sample_2"])
            self.assertEqual(config["abundance_weighted_samples"], ["sample_2", "sample_3"])
            self.assertEqual(config["reads_1"], {os.path.splitext(os.path.splitext(os.path.basename(s))[0])[0]: s for s in SAMPLE_READS_FORWARD.split(" ")})
            self.assertEqual(config["reads_2"], {os.path.splitext(os.path.splitext(os.path.basename(s))[0])[0]: s for s in SAMPLE_READS_REVERSE.split(" ")})
            self.assertEqual(config["assembly_strategy"], "megahit")

    def test_coassemble_singlem_inputs_files_of_paths(self):
        with in_tempdir():
            write_string_to_file(SAMPLE_READS_FORWARD, "sample_reads_forward")
            write_string_to_file(SAMPLE_READS_REVERSE, "sample_reads_reverse")
            write_string_to_file(GENOMES, "genomes")
            write_string_to_file(GENOME_TRANSCRIPTS, "genome_transcripts")
            write_string_to_file(SAMPLE_SINGLEM, "sample_singlem")

            cmd = (
                f"binchicken coassemble "
                f"--forward-list sample_reads_forward "
                f"--reverse-list sample_reads_reverse "
                f"--genomes-list genomes "
                f"--genome-transcripts-list genome_transcripts "
                f"--sample-singlem-list sample_singlem "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_singlem_inputs_dir(self):
        with in_tempdir():
            for sample in SAMPLE_READS_FORWARD.split(" ") + SAMPLE_READS_REVERSE.split(" "):
                extern.run(f"touch -t 200001011200 {sample}")

            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--sample-singlem-dir {SAMPLE_QUERY_DIR} "
                f"--sample-query-dir {SAMPLE_QUERY_DIR} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise " not in output)
            self.assertTrue("query_processing" in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_preclustered_dryrun(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--kmer-precluster always "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" not in output)
            self.assertTrue("sketch_samples" in output)
            self.assertTrue("distance_samples" in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" not in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_weighted_dryrun(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--abundance-weighted "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("abundance_weighting" in output)
            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("target_weighting" in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("qc_reads" not in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_metapackage_env_variable(self):
        with in_tempdir():
            os.environ['SINGLEM_METAPACKAGE_PATH'] = METAPACKAGE
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"singlem_appraise_filtered\" "
            )
            _ = subprocess.run(cmd, shell=True, check=True, capture_output=True, env=os.environ)

            appraise_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(appraise_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1.1", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "2", "3.28", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2.1", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3.1", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
                    ""
                ]
            )
            with open(appraise_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_remove_EIF(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--sample-singlem {SAMPLE_SINGLEM_EIF} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM_EIF} "
                f"--singlem-metapackage {METAPACKAGE_EIF} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"target_elusive\" "
            )
            output = extern.run(cmd)

            unbinned_path = os.path.join("test", "coassemble", "appraise", "unbinned.otu_table.tsv")
            self.assertTrue(os.path.exists(unbinned_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1.1", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "3", "4.92", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1.1", "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1.1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA", "1", "1.64", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1.1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis3", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2.1", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "4", "6.57", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2.1", "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC", "3", "4.92", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3.1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA", "6", "9.85", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3.1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis3", ""]),
                    ""
                ]
            )
            with open(unbinned_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_preclustered(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "4832"]),
                    ",".join(["sample_2", "3926"]),
                    ",".join(["sample_3", "3624"]),
                    ",".join(["sample_5", "3624"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
                self.assertEqual(expected, f.read())

            sketch_path = os.path.join("test", "coassemble", "sketch", "samples.sig")
            self.assertTrue(os.path.exists(sketch_path))

            distance_path = os.path.join("test", "coassemble", "sketch", "samples.csv")
            self.assertTrue(os.path.exists(distance_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "0,1"],
                    ["match", 2, "sample_1,sample_5", "0"],
                    ["match", 2, "sample_2,sample_5", "0"],
                    ["match", 2, "sample_3,sample_5", "3,4"],
                ],
                schema = ["style", "cluster_size", "samples", "target_ids"],
                orient="row",
            )
            observed = pl.read_csv(elusive_edges_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False)

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_3,sample_5",
                        "2",
                        "2",
                        "7248",
                        "sample_3,sample_5",
                        "coassembly_0"
                    ]),
                    "\t".join([
                        "sample_1,sample_2",
                        "2",
                        "2",
                        "8758",
                        "sample_1,sample_2,sample_5",
                        "coassembly_1"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_preclustered_target_taxa(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--taxa-of-interest p__Abyssobacteria "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--kmer-precluster always "
                f"--precluster-size 2 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sketch_path = os.path.join("test", "coassemble", "sketch", "samples.sig")
            self.assertTrue(os.path.exists(sketch_path))

            distance_path = os.path.join("test", "coassemble", "sketch", "samples.csv")
            self.assertTrue(os.path.exists(distance_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = "\n".join(
                [
                    "\t".join(["style", "cluster_size", "samples", "target_ids"]),
                    "\t".join(["match", "2", "sample_2,sample_5", "0"]),
                    ""
                ]
            )
            with open(elusive_edges_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_preclustered_single_assembly(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--single-assembly "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--max-recovery-samples 2 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "0,1"],
                    ["match", 2, "sample_1,sample_5", "0"],
                    ["match", 2, "sample_2,sample_5", "0"],
                    ["match", 2, "sample_3,sample_5", "3,4"],
                ],
                schema = ["style", "cluster_size", "samples", "target_ids"],
                orient="row",
            )
            observed = pl.read_csv(elusive_edges_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False)

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_5",
                        "1",
                        "3",
                        "3624",
                        "sample_3,sample_5",
                        "coassembly_0"
                    ]),
                    "\t".join([
                        "sample_3",
                        "1",
                        "2",
                        "3624",
                        "sample_3,sample_5",
                        "coassembly_1"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "2",
                        "3926",
                        "sample_1,sample_2",
                        "coassembly_2"
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "2",
                        "4832",
                        "sample_1,sample_2",
                        "coassembly_3"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_weighted(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--abundance-weighted "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_WEIGHT} "
                f"--genomes {GENOMES} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "4832"]),
                    ",".join(["sample_2", "3926"]),
                    ",".join(["sample_3", "3624"]),
                    ",".join(["sample_5", "3624"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
                self.assertEqual(expected, f.read())

            weights_path = os.path.join("test", "coassemble", "appraise", "weighted.otu_table.tsv")
            self.assertTrue(os.path.exists(weights_path))
            expected = pl.DataFrame([
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.126051],
                    ["S3.7.ribosomal_protein_S7", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 0.333333],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.479846],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.893421],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.022411],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.023992],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 0.169117],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.006483],
                ],
                schema = ["gene", "sequence", "weight"],
                orient="row",
            )
            observed = pl.read_csv(weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            target_weights_path = os.path.join("test", "coassemble", "target", "targets_weighted.tsv")
            self.assertTrue(os.path.exists(target_weights_path))
            expected = pl.DataFrame([
                    [0, 0.006483],
                    [1, 0.023992],
                    [2, 0.893421],
                    [3, 0.479846],
                    [4, 0.022411],
                    [5, 0.126051],
                    [6, 0.333333],
                    [7, 0.169117],
                ],
                schema = ["target", "weight"],
                orient="row",
            )
            observed = pl.read_csv(target_weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "1,3"],
                    ["match", 2, "sample_1,sample_5", "0,2"],
                    ["match", 2, "sample_2,sample_3", "5"],
                    ["match", 2, "sample_2,sample_5", "4,5"],
                    ["match", 2, "sample_3,sample_5", "5,7"],
                ],
                schema = ["style", "cluster_size", "samples", "target_ids"],
                orient="row",
            )
            observed = pl.read_csv(elusive_edges_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False)

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = pl.DataFrame([
                    ["sample_1,sample_5", 2, 0.899904, 8456, "sample_1,sample_5", "coassembly_0"],
                    ["sample_2,sample_3", 2, 0.126051, 7550, "sample_2,sample_3,sample_5", "coassembly_1"],
                ],
                schema = ["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"],
                orient="row",
            )
            observed = pl.read_csv(cluster_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

    def test_coassemble_weighted_specific_samples(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--abundance-weighted "
                f"--abundance-weighted-samples sample_1 sample_2 "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_WEIGHT} "
                f"--genomes {GENOMES} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))

            weights_path = os.path.join("test", "coassemble", "appraise", "weighted.otu_table.tsv")
            self.assertTrue(os.path.exists(weights_path))
            expected = pl.DataFrame([
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.008064],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.806451],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.479846],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.023992],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.039920],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.039920],
                ],
                schema = ["gene", "sequence", "weight"],
                orient="row",
            )
            observed = pl.read_csv(weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

    def test_coassemble_weighted_target_taxa(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--abundance-weighted "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_WEIGHT} "
                f"--taxa-of-interest p__Abyssobacteria "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            weights_path = os.path.join("test", "coassemble", "appraise", "weighted.otu_table.tsv")
            self.assertTrue(os.path.exists(weights_path))
            expected = pl.DataFrame([
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.022411],
                    ["S3.7.ribosomal_protein_S7", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", 0.048982],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 0.169117],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.006483],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.893421],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.126051],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.023992],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.479846],
                    ["S3.7.ribosomal_protein_S7", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 0.333333],

                ],
                schema = ["gene", "sequence", "weight"],
                orient="row",
            )
            observed = pl.read_csv(weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            target_weights_path = os.path.join("test", "coassemble", "target", "targets_weighted.tsv")
            self.assertTrue(os.path.exists(target_weights_path))
            expected = pl.DataFrame([
                    [0, 0.479846],
                    [1, 0.022411],
                    [2, 0.048982],
                    [3, 0.126051],
                    [4, 0.169117],
                ],
                schema = ["target", "weight"],
                orient="row",
            )
            observed = pl.read_csv(target_weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = "\n".join(
                [
                    "\t".join(["style", "cluster_size", "samples", "target_ids"]),
                    "\t".join(["match", "2", "sample_1,sample_2", "0"]),
                    "\t".join(["match", "2", "sample_2,sample_5", "1,3"]),
                    ""
                ]
            )
            with open(elusive_edges_path) as f:
                self.assertEqual(expected, f.read())

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = pl.DataFrame([
                    ["sample_1,sample_2", 2, 0.479846, 8758, "sample_1,sample_2", "coassembly_0"],
                ],
                schema = ["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"],
                orient="row",
            )
            observed = pl.read_csv(cluster_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

if __name__ == '__main__':
    unittest.main()
