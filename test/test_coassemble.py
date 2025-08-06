#!/usr/bin/env python3

import unittest
import os
import shutil
import gzip
from bird_tool_utils import in_tempdir
import extern
import subprocess
from snakemake.io import load_configfile
import polars as pl
from polars.testing import assert_frame_equal

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

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

PRECLUSTER_DISTANCES = os.path.join(path_to_data, "sketch", "sample_distances.csv")

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

    def test_coassemble_query_input_split(self):
        with in_tempdir():
            write_string_to_file(" ".join([f"sample_{i}.1.fq" for i in range(200)]), "sample_reads_forward")
            write_string_to_file(" ".join([f"sample_{i}.2.fq" for i in range(200)]), "sample_reads_reverse")
            write_string_to_file(GENOMES, "genomes")
            write_string_to_file(GENOME_TRANSCRIPTS, "genome_transcripts")
            with open("sample_size.csv", "w") as f:
                f.write("\n".join([",".join([f"sample_{i}", "100"]) for i in range(200)]))

            samples = [f"sample_{i}" for i in range(200)]
            query_dir = "sample_queries"
            os.makedirs(query_dir)
            singlem_dir = "sample_singlem"
            os.makedirs(singlem_dir)
            for sample in samples:
                with open(os.path.join(query_dir, sample + "_query.otu_table.tsv"), "w") as f:
                    f.write("\t".join(["query_name", "query_sequence", "divergence", "num_hits", "coverage", "sample", "marker", "hit_sequence", "taxonomy"]))
                    f.write("\n")
                    f.write("\t".join([
                        sample,
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "10",
                        "1",
                        "1.64",
                        "GB_GCA_013286235.1",
                        "S3.7.ribosomal_protein_S7",
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335"
                    ]))

                with open(os.path.join(singlem_dir, sample + "_read.otu_table.tsv"), "w") as f:
                    f.write("\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy"]))
                    f.write("\n")
                    f.write("\t".join([
                        "S3.7.ribosomal_protein_S7",
                        sample,
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "3",
                        "4.92",
                        "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335"
                    ]))

            cmd = (
                f"binchicken coassemble "
                f"--forward-list sample_reads_forward "
                f"--reverse-list sample_reads_reverse "
                f"--genomes-list genomes "
                f"--genome-transcripts-list genome_transcripts "
                f"--sample-query-dir {query_dir} "
                f"--sample-singlem-dir {singlem_dir} "
                f"--sample-read-size sample_size.csv "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--snakemake-args \"cluster_graph\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    ""
                ]
            )
            with open(binned_path) as f:
                self.assertEqual(expected, f.read())

            unbinned_path = os.path.join("test", "coassemble", "appraise", "unbinned.otu_table.tsv")
            self.assertTrue(os.path.exists(unbinned_path))
            expected = "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in",]) + "\n" + "\n".join(
                [
                    "\t".join([
                        "S3.7.ribosomal_protein_S7",
                        sample,
                        "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT",
                        "3",
                        "4.92",
                        "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335",
                        "",
                    ])
                    for sample in samples
                ]
            ) + "\n"
            with open(unbinned_path) as f:
                self.assertEqual(expected, f.read())

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

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
                        "sample_1"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "3",
                        "3926",
                        "sample_1,sample_2,sample_3",
                        "sample_2"
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
            write_string_to_file("sample_1\nsample_3", "anchor_samples")

            cmd = (
                f"binchicken coassemble "
                f"--forward-list sample_reads_forward "
                f"--reverse-list sample_reads_reverse "
                f"--anchor-samples-list anchor_samples "
                f"--genomes-list genomes "
                f"--genome-transcripts-list genome_transcripts "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--coassembly-samples-list coassembly_samples "
                f"--exclude-coassemblies-list exclude_coassemblies "
                f"--abundance-weighted-samples-list abundance_samples "
                f"--assemble-unmapped "
                f"--assembly-strategy megahit "
                f"--output test "
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
            self.assertEqual(config["anchor_samples"], ["sample_1", "sample_3"])

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

    def test_coassemble_anchored(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--anchor-samples sample_3 "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--prodigal-meta "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "coassemble", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))

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
                f"--snakemake-args \"singlem_appraise_filtered\" "
            )
            _ = subprocess.run(cmd, shell=True, check=True, capture_output=True, env=os.environ)

            appraise_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(appraise_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "2", "3.28", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3", "TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC", "1", "1.64", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235", "GB_GCA_013286235.1_protein"]),
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
                f"--snakemake-args \"target_elusive\" "
            )
            output = extern.run(cmd)

            unbinned_path = os.path.join("test", "coassemble", "appraise", "unbinned.otu_table.tsv")
            self.assertTrue(os.path.exists(unbinned_path))
            expected = "\n".join(
                [
                    "\t".join(["gene", "sample", "sequence", "num_hits", "coverage", "taxonomy", "found_in"]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "3", "4.92", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA", "1", "1.64", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_1", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis3", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "4", "6.57", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter sp018688335", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_2", "TACCAGGTCCCGGTCGAGGTCCGTCCGATCCGCCAGACGACGCTCGCCCTGCGCTGGCTC", "3", "4.92", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATA", "6", "9.85", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis2", ""]),
                    "\t".join(["S3.7.ribosomal_protein_S7", "sample_3", "TATCAGGTGCCTATTGAGGTAAGACCTGAAAGAAGACAGACTTTAGCGCTTCGCTGGATC", "5", "8.21", "Root; d__Bacteria; p__Actinobacteriota; c__Actinomycetia; o__Mycobacteriales; f__Mycobacteriaceae; g__Nocardia; s__Nocardia grenadensis3", ""]),
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
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_2,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_3,sample_5", "1779956245962588283,5034568815038442683"],
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
                    "\t".join(["match", "2", "sample_2,sample_5", "2416810233479675551"]),
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
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_2,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_3,sample_5", "1779956245962588283,5034568815038442683"],
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
                        "sample_5"
                    ]),
                    "\t".join([
                        "sample_3",
                        "1",
                        "2",
                        "3624",
                        "sample_3,sample_5",
                        "sample_3"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "2",
                        "3926",
                        "sample_1,sample_2",
                        "sample_2"
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "2",
                        "4832",
                        "sample_1,sample_2",
                        "sample_1"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_preclustered_distances_provided(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--precluster-distances {PRECLUSTER_DISTANCES} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)

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
            self.assertFalse(os.path.exists(sketch_path))

            distance_path = os.path.join("test", "coassemble", "sketch", "samples.csv")
            self.assertTrue(os.path.exists(distance_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_2,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_3,sample_5", "1779956245962588283,5034568815038442683"],
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

    def test_coassemble_preclustered_anchored(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--anchor-samples sample_1 "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--output test "
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
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
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
                        "sample_1,sample_2",
                        "2",
                        "2",
                        "8758",
                        "sample_1,sample_2,sample_5",
                        "coassembly_0"
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
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.094538],
                    ["S3.7.ribosomal_protein_S7", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 0.083333],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.239923],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.446710],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.011205],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.011996],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 0.084558],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.003241],
                ],
                schema = ["gene", "sequence", "weight"],
                orient="row",
            )
            observed = pl.read_csv(weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            target_weights_path = os.path.join("test", "coassemble", "target", "targets_weighted.tsv")
            self.assertTrue(os.path.exists(target_weights_path))
            expected = pl.DataFrame([
                    [0, 0.003241],
                    [1, 0.011996],
                    [2, 0.446710],
                    [3, 0.239923],
                    [4, 0.011205],
                    [5, 0.094538],
                    [6, 0.083333],
                    [7, 0.084558],
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
                    ["sample_1,sample_5", 2, 0.003241+0.446710, 8456, "sample_1,sample_5", "coassembly_0"],
                    ["sample_2,sample_3", 2, 0.094538, 7550, "sample_2,sample_3,sample_5", "coassembly_1"],
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
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))

            weights_path = os.path.join("test", "coassemble", "appraise", "weighted.otu_table.tsv")
            self.assertTrue(os.path.exists(weights_path))
            expected = pl.DataFrame([
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.004032],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.403225],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.479846],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.023992],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.019960],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.019960],
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
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            weights_path = os.path.join("test", "coassemble", "appraise", "weighted.otu_table.tsv")
            self.assertTrue(os.path.exists(weights_path))
            expected = pl.DataFrame([
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCG---------------GTTATATATCGAAAGCTCATGCGGCCATATCG", 0.094538],
                    ["S3.7.ribosomal_protein_S7", "ATCGACTGACTTGATCGATCTTTGACGACGAGAGAGAGAGCGACGCGCCGAGAGGTTTCA", 0.083333],
                    ["S3.7.ribosomal_protein_S7", "AAAAGTCCTGATCGTAGCTAATAATATTATGCGTACGTCAGTACGTACTGACTGACGTAA", 0.239923],
                    ["S3.7.ribosomal_protein_S7", "AGCGTCGAGCGATCGATCGTACGTAGCGGGGATCGTATTTACTATCTACTAACGAGCAAA", 0.446710],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAGTTAGGAAAGCCCCCGGAGTTAGCTA", 0.011205],
                    ["S3.7.ribosomal_protein_S7", "TGACTAGCTGGGCTAGCTATATTCTTTTTACGAGCGCGAGGAAAGCGACAGCGGCCAGGC", 0.011996],
                    ["S3.7.ribosomal_protein_S7", "TACGAGCGGATCGTGCACGTAGTCAGTCGTTATATATCGAAAGCTCATGCGGCCATATCG", 0.084558],
                    ["S3.7.ribosomal_protein_S7", "ATGACTAGTCATAGCTAGATTTGAGGCAGCAGGAAGCAGCATAGCTACGACACACCCCCC", 0.003241],
                    ["S3.7.ribosomal_protein_S7", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", 0.024491],
                ],
                schema = ["gene", "sequence", "weight"],
                orient="row",
            )
            observed = pl.read_csv(weights_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

            target_weights_path = os.path.join("test", "coassemble", "target", "targets_weighted.tsv")
            self.assertTrue(os.path.exists(target_weights_path))
            expected = pl.DataFrame([
                    [0, 0.239923],
                    [1, 0.011205],
                    [2, 0.024491],
                    [3, 0.094538],
                    [4, 0.084558],
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
                    ["sample_1,sample_2", 2, 0.239923, 8758, "sample_1,sample_2", "coassembly_0"],
                ],
                schema = ["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"],
                orient="row",
            )
            observed = pl.read_csv(cluster_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False, atol=1e-5)

    def test_coassemble_sra_download(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output_comb = extern.run(cmd)

            output_sra = output_comb.split("Building DAG of jobs...")[1]
            self.assertTrue("download_sra" in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output = output_comb.split("Building DAG of jobs...")[2]
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("genome_transcripts" in output)
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
            self.assertTrue("qc_reads" in output)
            self.assertTrue("collect_genomes" not in output)
            self.assertTrue("map_reads" not in output)
            self.assertTrue("finish_mapping" not in output)
            self.assertTrue("aviary_commands" in output)

    def test_coassemble_sra_download_mock(self):
        with in_tempdir():
            cmd = (
                f"binchicken coassemble "
                f"--forward SRR3309137 SRR8334323 SRR8334324 "
                f"--sra "
                f"--output test "
                f"--snakemake-args \" --config mock_sra=True\" "
            )
            extern.run(cmd)

            sra_f0_path = os.path.join("test", "coassemble", "sra", "SRR3309137_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f0_path))
            with gzip.open(sra_f0_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" in file)
                self.assertTrue("@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" not in file)

            sra_r0_path = os.path.join("test", "coassemble", "sra", "SRR3309137_2.fastq.gz")
            self.assertTrue(os.path.exists(sra_r0_path))
            with gzip.open(sra_r0_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" in file)
                self.assertTrue("@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" not in file)

            sra_f1_path = os.path.join("test", "coassemble", "sra", "SRR8334323_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f1_path))
            with gzip.open(sra_f1_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SEQ_ID.1" in file)

            sra_f2_path = os.path.join("test", "coassemble", "qc", "SRR8334323_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f2_path))
            with gzip.open(sra_f2_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SEQ_ID.1" not in file)

            recover_path = os.path.join("test", "coassemble", "target", "unused_samples.tsv")
            self.assertTrue(os.path.exists(recover_path))
            with open(recover_path) as f:
                file = [l.strip() for l in f.readlines()]
                self.assertTrue("SRR3309137" in file)
                self.assertTrue("SRR8334323" in file)
                self.assertTrue("SRR8334324" in file)

    def test_coassemble_hierarchy(self):
        with in_tempdir():
            local_samples_forward = [
                os.path.abspath("sample_1.1.fq"),
                os.path.abspath("different_2.1.fq"),
                os.path.abspath("simple_3.1.fq"),
                os.path.abspath("sam_4.1.fq"),
            ]
            for old, new in zip(SAMPLE_READS_FORWARD_EMPTY.split(" "), local_samples_forward):
                shutil.copy(old, new)

            local_samples_reverse = [
                os.path.abspath("sample_1.2.fq"),
                os.path.abspath("different_2.2.fq"),
                os.path.abspath("simple_3.2.fq"),
                os.path.abspath("sam_4.2.fq"),
            ]
            for old, new in zip(SAMPLE_READS_REVERSE_EMPTY.split(" "), local_samples_reverse):
                shutil.copy(old, new)

            cmd = (
                f"binchicken coassemble "
                f"--forward {' '.join(local_samples_forward)} "
                f"--reverse {' '.join(local_samples_reverse)} "
                f"--genomes {TWO_GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--assemble-unmapped "
                f"--unmapping-max-identity 99 "
                f"--unmapping-max-alignment 90 "
                f"--prodigal-meta "
                f"--run-qc "
                f"--file-hierarchy always "
                f"--file-hierarchy-depth 3 "
                f"--file-hierarchy-chars 3 "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            for sample in local_samples_forward:
                sample_name = os.path.basename(sample).split(".")[0]
                if len(sample_name) > 6:
                    subpath = os.path.join(sample_name[0:3], sample_name[3:6], sample_name[6:9])
                else:
                    subpath = os.path.join(sample_name[0:3], sample_name[3:6])

                pipe_path = os.path.join("test", "coassemble", "pipe", subpath, sample_name + "_read.otu_table.tsv")
                self.assertTrue(os.path.exists(pipe_path))

                qc_for_path = os.path.join("test", "coassemble", "qc", subpath, sample_name + "_1.fastq.gz")
                self.assertTrue(os.path.exists(qc_for_path))

                qc_rev_path = os.path.join("test", "coassemble", "qc", subpath, sample_name + "_2.fastq.gz")
                self.assertTrue(os.path.exists(qc_rev_path))

                map_for_path = os.path.join("test", "coassemble", "mapping", subpath, sample_name + "_unmapped.1.fq.gz")
                self.assertTrue(os.path.exists(map_for_path))

                map_rev_path = os.path.join("test", "coassemble", "mapping", subpath, sample_name + "_unmapped.2.fq.gz")
                self.assertTrue(os.path.exists(map_rev_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "4832"]),
                    ",".join(["different_2", "3926"]),
                    ",".join(["simple_3", "3624"]),
                    ",".join(["sam_4", "604"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
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
                        "different_2,sample_1",
                        "2",
                        "3",
                        "8758",
                        "different_2,sample_1,simple_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_sample_naming(self):
        with in_tempdir():
            local_samples_forward = [
                os.path.abspath("sample_1.1.fq"),
                os.path.abspath("sample_2_R1.fq"),
                os.path.abspath("sample_3_1.fq"),
                os.path.abspath("sample_R1_4.fq"),
            ]
            for old, new in zip(SAMPLE_READS_FORWARD_EMPTY.split(" "), local_samples_forward):
                shutil.copy(old, new)

            local_samples_reverse = [
                os.path.abspath("sample_1.2.fq"),
                os.path.abspath("sample_2_R2.fq"),
                os.path.abspath("sample_3_2.fq"),
                os.path.abspath("sample_R2_4.fq"),
            ]
            for old, new in zip(SAMPLE_READS_REVERSE_EMPTY.split(" "), local_samples_reverse):
                shutil.copy(old, new)

            cmd = (
                f"binchicken coassemble "
                f"--forward {' '.join(local_samples_forward)} "
                f"--reverse {' '.join(local_samples_reverse)} "
                f"--genomes {TWO_GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--prodigal-meta "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sample_1_path = os.path.join("test", "coassemble", "pipe", "sample_1_read.otu_table.tsv")
            self.assertTrue(os.path.exists(sample_1_path))

            sample_2_path = os.path.join("test", "coassemble", "pipe", "sample_2_read.otu_table.tsv")
            self.assertTrue(os.path.exists(sample_2_path))

            sample_3_path = os.path.join("test", "coassemble", "pipe", "sample_3_read.otu_table.tsv")
            self.assertTrue(os.path.exists(sample_3_path))

            sample_4_path = os.path.join("test", "coassemble", "pipe", "sample_read.otu_table.tsv")
            self.assertTrue(os.path.exists(sample_4_path))

            read_size_path = os.path.join("test", "coassemble", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "4832"]),
                    ",".join(["sample_2", "3926"]),
                    ",".join(["sample_3", "3624"]),
                    ",".join(["sample", "604"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
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

if __name__ == '__main__':
    unittest.main()
