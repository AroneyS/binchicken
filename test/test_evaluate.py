#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
from snakemake.io import load_configfile

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_UNBINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned.otu_table.tsv")
MOCK_BINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "binned.otu_table.tsv")
MOCK_TARGETS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "targets.tsv")
MOCK_ELUSIVE_EDGES = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_edges.tsv")
MOCK_ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_clusters.tsv")
MOCK_SUMMARY = os.path.join(MOCK_COASSEMBLE, "coassemble", "summary.tsv")

MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])
MOCK_GENOMES = " ".join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_1.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_2.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_3.fna"),
])

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

class Tests(unittest.TestCase):
    def test_evaluate(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--checkm-version 2 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "statistic",
                        "within",
                        "match",
                        "nonmatch",
                        "total",
                        "match_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "recovery",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_bin_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_unbin_sequences",
                        "recovery",
                        "0",
                        "3",
                        "3",
                        "0.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "novel_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-0_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.18", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_specified_files(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("prodigal_bins" in output)
            self.assertTrue("singlem_pipe_bins" in output)
            self.assertTrue("singlem_summarise_bins" in output)
            self.assertTrue("cluster_original_bins" not in output)
            self.assertTrue("cluster_updated_bins" not in output)
            self.assertTrue("summarise_clusters" not in output)
            self.assertTrue("evaluate" in output)
            self.assertTrue("evaluate_plots" in output)

    def test_evaluate_genome_input(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--new-genomes {MOCK_GENOMES} "
                f"--prodigal-meta "
                f"--coassembly-run coassembly_0 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "statistic",
                        "within",
                        "match",
                        "nonmatch",
                        "total",
                        "match_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_bin_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_unbin_sequences",
                        "recovery",
                        "0",
                        "3",
                        "3",
                        "0.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "novel_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_1_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_3_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_3_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.14", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_genome_files_input(self):
        with in_tempdir():
            write_string_to_file(MOCK_GENOMES, "genomes")

            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--new-genomes-list genomes "
                f"--prodigal-meta "
                f"--coassembly-run coassembly_0 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("prodigal_bins" in output)
            self.assertTrue("singlem_pipe_bins" in output)
            self.assertTrue("singlem_summarise_bins" in output)
            self.assertTrue("cluster_original_bins" not in output)
            self.assertTrue("cluster_updated_bins" not in output)
            self.assertTrue("summarise_clusters" not in output)
            self.assertTrue("evaluate" in output)
            self.assertTrue("evaluate_plots" in output)

    def test_evaluate_default_config(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["checkm_version"], 2)
            self.assertEqual(config["min_completeness"], 70)
            self.assertEqual(config["max_contamination"], 10)

    def test_evaluate_metapackage_env_variable(self):
        with in_tempdir():
            os.environ['SINGLEM_METAPACKAGE_PATH'] = METAPACKAGE
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \"singlem_summarise_bins\" "
            )
            import subprocess
            _ = subprocess.run(cmd, shell=True, check=True, capture_output=True, env=os.environ)

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-0_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.18", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_cluster(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--cluster "
                f"--genomes {TWO_GENOMES} "
                f"--checkm-version 2 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))

            cluster_path = os.path.join("test", "evaluate", "evaluate", "cluster_stats.csv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    ",".join([
                        "original",
                        "2",
                    ]),
                    ",".join([
                        "coassembly_0",
                        "4",
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())


if __name__ == '__main__':
    unittest.main()
