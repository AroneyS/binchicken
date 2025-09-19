#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
import subprocess
from snakemake.io import load_configfile
import re
import polars as pl

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
SAMPLE_READS_FORWARD_NO_TWO = " ".join([
    os.path.join(path_to_data, "sample_1.1.fq"),
    os.path.join(path_to_data, "sample_3.1.fq"),
])
SAMPLE_READS_REVERSE_NO_TWO = " ".join([
    os.path.join(path_to_data, "sample_1.2.fq"),
    os.path.join(path_to_data, "sample_3.2.fq"),
])

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_COASSEMBLE_MINIMAL = os.path.join(path_to_data, "mock_coassemble_minimal")
MOCK_UNBINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned.otu_table.tsv")
MOCK_UNBINNED_BIASED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned_biased.otu_table.tsv")
MOCK_BINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "binned.otu_table.tsv")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])
MOCK_GENOMES = " ".join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_1.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_2.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_3.fna"),
])
MOCK_GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "coassemble", "summarise", "bins_summarised_mock.otu_table.tsv")
ELUSIVE_CLUSTERS = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_clusters_iterate.tsv")])

SAMPLE_READ_SIZE = os.path.join(MOCK_COASSEMBLE, "coassemble", "read_size2.csv")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_3_read.otu_table.tsv"),
    ])
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "coassemble", "summarise", "bins_summarised.otu_table2.tsv")

BIN_INFO_COLUMNS = {
    "Bin Id": str,
    "Marker lineage": str,
    "# genomes": str,
    "# markers": str,
    "# marker sets": str,
    "Completeness (CheckM1)": str,
    "Contamination (CheckM1)": str,
    "Strain heterogeneity": str,
    "Genome size (bp)": str,
    "# ambiguous bases": str,
    "# scaffolds": str,
    "# contigs": str,
    "N50 (scaffolds)": str,
    "N50 (contigs)": str,
    "Mean scaffold length (bp)": str,
    "Mean contig length (bp)": str,
    "Longest scaffold (bp)": str,
    "Longest contig (bp)": str,
    "GC": str,
    "GC std (scaffolds > 1kbp)": str,
    "Coding density": str,
    "Translation table": str,
    "# predicted genes": str,
    "0": str,
    "1": str,
    "2": str,
    "3": str,
    "4": str,
    "5+": str,
    "Completeness (CheckM2)": str,
    "Contamination (CheckM2)": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

ELUSIVE_CLUSTERS_COLUMNS = {
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
}

CHECKM2_QUALITY_COLUMNS = {
    "Name": str,
    "Completeness": str,
    "Contamination": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

class Tests(unittest.TestCase):
    def test_iterate(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--iteration 0 "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            extern.run(cmd)

            new_bin_0_path = os.path.join("test", "recovered_bins", "iteration_0-coassembly_0-0.fna")
            self.assertTrue(os.path.exists(new_bin_0_path))
            with open(new_bin_0_path) as f:
                file = f.read()
                self.assertTrue(">k141_1363016" in file)
                self.assertTrue(">k141_177318" not in file)

            new_bin_1_path = os.path.join("test", "recovered_bins", "iteration_0-coassembly_0-1.fna")
            self.assertTrue(os.path.exists(new_bin_1_path))

            new_bin_2_path = os.path.join("test", "recovered_bins", "iteration_0-coassembly_0-2.fna")
            self.assertFalse(os.path.exists(new_bin_2_path))

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            original_bin_singlem_path = os.path.join("test", "coassemble", "pipe", os.path.splitext(os.path.basename(GENOMES.split(" ")[0]))[0]+"_bin.otu_table.tsv")
            self.assertFalse(os.path.exists(original_bin_singlem_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            with open(binned_path) as f:
                file = f.read()
                self.assertTrue("TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG" in file)
                self.assertTrue("GB_GCA_013286235.1_protein" in file)
                self.assertTrue("TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT" in file)
                self.assertTrue("iteration_0-coassembly_0-0_protein" in file)
                # Check that header is not repeated
                header = "gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\tfound_in"
                self.assertFalse(header in file.replace(header, "", 1))

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
                        "1",
                        "8456",
                        "sample_1,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_iterate_minimal(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            self.assertTrue("count_bp_reads" not in output)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            NEW_GENOMES = " ".join([
                os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "iteration_0-coassembly_0-0.fna"),
                os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "iteration_0-coassembly_0-1.fna"),
            ])
            genomes = {
                os.path.splitext(os.path.basename(g))[0]: g.replace(MOCK_COASSEMBLE + "/coassemble/coassembly_0/recover/bins/final_bins/", os.getcwd() + "/test/recovered_bins/")
                for g in (GENOMES + " " + NEW_GENOMES).split(" ")
                }
            self.assertEqual(genomes, config["genomes"])

            reads_1 = {
                os.path.splitext(os.path.basename(r))[0].removesuffix(".1"): r
                for r in SAMPLE_READS_FORWARD.split(" ")
                }
            self.assertEqual(reads_1, config["reads_1"])

            reads_2 = {
                os.path.splitext(os.path.basename(r))[0].removesuffix(".2"): r
                for r in SAMPLE_READS_REVERSE.split(" ")
                }
            self.assertEqual(reads_2, config["reads_2"])

            exclude_coassemblies = ["sample_0,sample_1", "sample_1,sample_2"]
            self.assertEqual(exclude_coassemblies, config["exclude_coassemblies"])

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
                        "1",
                        "8456",
                        "sample_1,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

            cumulative_coassemblies_path = os.path.join("test", "coassemble", "target", "cumulative_coassemblies.tsv")
            self.assertTrue(os.path.exists(cumulative_coassemblies_path))
            expected = "\n".join(["sample_0,sample_1", "sample_1,sample_2", ""])
            with open(cumulative_coassemblies_path) as f:
                self.assertEqual(expected, f.read())

    def test_iterate_minimal_input(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--coassemble-output {MOCK_COASSEMBLE_MINIMAL} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            self.assertTrue("count_bp_reads" not in output)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            NEW_GENOMES = [
                os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "iteration_0-coassembly_0-0.fna"),
                os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "iteration_0-coassembly_0-1.fna"),
            ]
            genomes = {
                os.path.splitext(os.path.basename(g))[0]: g.replace(MOCK_COASSEMBLE + "/coassemble/coassembly_0/recover/bins/final_bins/", os.getcwd() + "/test/recovered_bins/")
                for g in NEW_GENOMES
                }
            self.assertEqual(genomes, config["genomes"])

            reads_1 = {
                os.path.splitext(os.path.basename(r))[0].removesuffix(".1"): r
                for r in SAMPLE_READS_FORWARD.split(" ")
                }
            self.assertEqual(reads_1, config["reads_1"])

            reads_2 = {
                os.path.splitext(os.path.basename(r))[0].removesuffix(".2"): r
                for r in SAMPLE_READS_REVERSE.split(" ")
                }
            self.assertEqual(reads_2, config["reads_2"])

            exclude_coassemblies = ["sample_1,sample_2"]
            self.assertEqual(exclude_coassemblies, config["exclude_coassemblies"])

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
                        "1",
                        "8456",
                        "sample_1,sample_3",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

            cumulative_coassemblies_path = os.path.join("test", "coassemble", "target", "cumulative_coassemblies.tsv")
            self.assertTrue(os.path.exists(cumulative_coassemblies_path))
            expected = "\n".join(["sample_1,sample_2", ""])
            with open(cumulative_coassemblies_path) as f:
                self.assertEqual(expected, f.read())

    def test_iterate_genome_input(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--new-genomes {MOCK_GENOMES} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--prodigal-meta "
            )
            extern.run(cmd)

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))
            expected = "\n".join(["\t".join([g, os.path.basename(g)]) for g in MOCK_GENOMES.split(" ")])
            with open(bin_provenance_path) as f:
                self.assertEqual(expected, f.read())

            original_bin_singlem_path = os.path.join("test", "coassemble", "pipe", os.path.splitext(os.path.basename(GENOMES.split(" ")[0]))[0]+"_bin.otu_table.tsv")
            self.assertTrue(os.path.exists(original_bin_singlem_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            with open(binned_path) as f:
                file = f.read()
                self.assertTrue("TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC" in file)
                self.assertTrue("GB_GCA_013286235.1_protein" in file)
                self.assertTrue("TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT" in file)
                self.assertTrue("bin_1_protein" in file)

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

    def test_iterate_genome_files_input(self):
        with in_tempdir():
            write_string_to_file(MOCK_GENOMES, "genomes")
            cmd = (
                f"binchicken iterate "
                f"--new-genomes-list genomes "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" in output)
            self.assertTrue("singlem_pipe_genomes" in output)
            self.assertTrue("singlem_summarise_genomes" in output)
            self.assertTrue("update_appraise" in output)
            self.assertTrue("singlem_appraise_filtered" in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("cluster_graph" in output)
            self.assertTrue("aviary_commands" in output)
            self.assertTrue("summary" in output)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            genomes_config = {
                os.path.splitext(os.path.basename(g))[0]: g
                for g in (GENOMES + " " + MOCK_GENOMES).split(" ")
                }
            self.assertEqual(genomes_config, config["genomes"])

            binned_prior_path = os.path.join("test", "coassemble", "appraise", "binned_prior.otu_table.tsv")
            self.assertTrue(os.path.islink(binned_prior_path))
            self.assertEqual(MOCK_BINNED, os.path.realpath(binned_prior_path))

            unbinned_prior_path = os.path.join("test", "coassemble", "appraise", "unbinned_prior.otu_table.tsv")
            self.assertTrue(os.path.islink(unbinned_prior_path))
            self.assertEqual(MOCK_UNBINNED, os.path.realpath(unbinned_prior_path))

    def test_iterate_default_config(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--exclude-coassemblies sample_4,sample_5 "
                f"--output test "
                f"--dryrun "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["taxa_of_interest"], "")
            self.assertEqual(config["assemble_unmapped"], False)
            self.assertEqual(config["aviary_assemble_threads"], 64)
            self.assertEqual(config["aviary_assemble_memory"], 500)
            self.assertEqual(config["aviary_recover_threads"], 32)
            self.assertEqual(config["aviary_recover_memory"], 250)
            self.assertEqual(config["exclude_coassemblies"], ["sample_0,sample_1", "sample_1,sample_2", "sample_4,sample_5"])

            self.assertTrue("Evaluating bins using CheckM2 with completeness >= 70 and contamination <= 10" in output)

    def test_iterate_genome_singlem(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--iteration 0 "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--exclude-coassemblies sample_1,sample_2 "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            self.assertTrue(re.search(r"singlem_pipe_genomes\s+2", output))

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

    def test_iterate_new_genome_singlem(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--new-genomes {MOCK_GENOMES} "
                f"--new-genome-singlem {MOCK_GENOME_SINGLEM} "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            self.assertFalse(re.search(r"singlem_pipe_genomes\s", output))

            original_bin_singlem_path = os.path.join("test", "coassemble", "pipe", os.path.splitext(os.path.basename(GENOMES.split(" ")[0]))[0]+"_bin.otu_table.tsv")
            self.assertFalse(os.path.exists(original_bin_singlem_path))

            new_bin_singlem_path = os.path.join("test", "coassemble", "pipe", os.path.splitext(os.path.basename(MOCK_GENOMES.split(" ")[0]))[0]+"_bin.otu_table.tsv")
            self.assertFalse(os.path.exists(new_bin_singlem_path))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            with open(binned_path) as f:
                file = f.read()
                self.assertTrue("TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG" in file)
                self.assertTrue("GB_GCA_013286235.1_protein" in file)
                self.assertTrue("TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT" in file)
                self.assertTrue("bin_1_protein" in file)

    def test_iterate_genome_singlem_and_new(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--new-genomes {MOCK_GENOMES} "
                f"--new-genome-singlem {MOCK_GENOME_SINGLEM} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            self.assertFalse(re.search(r"singlem_pipe_genomes\s", output))

            binned_path = os.path.join("test", "coassemble", "appraise", "binned.otu_table.tsv")
            self.assertTrue(os.path.exists(binned_path))
            with open(binned_path) as f:
                file = f.read()
                self.assertTrue("TTCCAGGTGCCTACCGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTC" in file)
                self.assertTrue("GB_GCA_013286235.1_protein" in file)
                self.assertTrue("TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT" in file)
                self.assertTrue("bin_1_protein" in file)

    def test_iterate_missing_samples(self):
        with in_tempdir():
            cmd = (
                f"binchicken iterate "
                f"--forward {SAMPLE_READS_FORWARD_NO_TWO} "
                f"--reverse {SAMPLE_READS_REVERSE_NO_TWO} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--coassemble-unbinned {MOCK_UNBINNED_BIASED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--genomes {GENOMES} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            extern.run(cmd)

            bin_provenance_path = os.path.join("test", "previous_bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            elusive_clusters_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(elusive_clusters_path))
            expected = "\n".join(
                [
                    "\t".join(["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"]),
                    "\t".join(["sample_1,sample_3", "2", "1", "2869", "sample_1,sample_3", "coassembly_0"]),
                    ""
                ]
            )
            with open(elusive_clusters_path) as f:
                self.assertEqual(expected, f.read())

    def test_iterate_prefers_recovered_bins(self):
        with in_tempdir():
            # elusive clusters
            elusive_clusters_path = os.path.join("prior_coassemble", "coassemble", "target", "elusive_clusters.tsv")
            os.makedirs(os.path.dirname(elusive_clusters_path), exist_ok=True)
            (
                pl.DataFrame([
                        ["sample_1,sample_2", "2", "2", "2869", "sample_1,sample_2", "coassembly_0"],
                        ["sample_1,sample_3", "2", "2", "2869", "sample_1,sample_3", "coassembly_1"],
                        ],
                    orient="row", schema=ELUSIVE_CLUSTERS_COLUMNS
                    )
                .write_csv(elusive_clusters_path, separator="\t")
            )

            unbinned_path = os.path.join("prior_coassemble", "coassemble", "appraise", "unbinned.otu_table.tsv")
            binned_path = os.path.join("prior_coassemble", "coassemble", "appraise", "binned.otu_table.tsv")
            os.makedirs(os.path.dirname(unbinned_path), exist_ok=True)
            read_size_path = os.path.join("prior_coassemble", "coassemble", "read_size.csv")

            with open(unbinned_path, "w") as f:
                pass
            with open(binned_path, "w") as f:
                pass
            with open(read_size_path, "w") as f:
                pass

            # recovered_bins outputs
            output_bins_path = os.path.join("prior_coassemble", "recovered")
            os.makedirs(os.path.join(output_bins_path, "bins"), exist_ok=True)

            genome1_path = os.path.join(output_bins_path, "bins", "binchicken_coassembly_0.1.fna")
            with open(genome1_path, "w") as f:
                f.write(">contigA\nACGT\n")

            genome2_path = os.path.join(output_bins_path, "bins", "binchicken_coassembly_0.2.fna")
            with open(genome2_path, "w") as f:
                f.write(">contigB\nTGCA\n")

            bin_info = pl.DataFrame([
                    [
                        "binchicken_coassembly_0.1",
                        "p__Bacteroidetes (UID2605)", "350", "314", "208", "86.44", "1.8", "12.5", "1705643", "0", "315", "315", "6824", "6824", "5414", "5414", "25260", "25260", "59.9", "3.41", "90.22", "11",
                        "1669", "37", "269", "8", "0", "0", "0",
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "binchicken_coassembly_0.2",
                        "o__Clostridiales (UID1212)", "172", "257", "149", "44.16", "1.03", "33.33", "980304", "0", "399", "399", "2442", "2442", "2456", "2456", "7735", "7735", "60.2", "2.97", "90.94", "11",
                        "1244", "143", "111", "3", "0", "0", "0",
                        "36.49", "0.24", "Neural Network (Specific Model)", "11", "0.909", "2442", "239.274115755627", "980304", "0.6", "1244", "399", "7735", ""
                    ],
                ], orient="row", schema=BIN_INFO_COLUMNS)
            (
                bin_info
                .write_csv(os.path.join(output_bins_path, "bin_info.tsv"), separator="\t")
            )

            (
                bin_info
                .rename({
                    "Bin Id": "Name",
                    "Completeness (CheckM2)": "Completeness",
                    "Contamination (CheckM2)": "Contamination",
                    })
                .select(CHECKM2_QUALITY_COLUMNS.keys())
                .write_csv(os.path.join(output_bins_path, "quality_report.tsv"), separator="\t")
            )

            config_path = os.path.join("prior_coassemble", "config.yaml")
            with open(config_path, "w") as f:
                f.write("\n".join([
                    "no_genomes: true",
                    f"reads_1: {{'sample_1': '{path_to_data}/sample_1.1.fq', 'sample_2': '{path_to_data}/sample_2.1.fq', 'sample_3': '{path_to_data}/sample_3.1.fq'}}",
                    f"reads_2: {{'sample_1': '{path_to_data}/sample_1.2.fq', 'sample_2': '{path_to_data}/sample_2.2.fq', 'sample_3': '{path_to_data}/sample_3.2.fq'}}",
                    "",
                ]))

            # competing coassembly_1 outputs
            source_bins_path = os.path.join("prior_coassemble", "coassemble", "coassemble", "coassembly_1", "recover", "bins")
            os.makedirs(os.path.join(source_bins_path, "final_bins"), exist_ok=True)

            source_genome1_path = os.path.join(source_bins_path, "final_bins", "metabat2_refined_bins.tsv.005.fna")
            with open(source_genome1_path, "w") as f:
                f.write(">contigA\nACGT\n")

            source_genome2_path = os.path.join(source_bins_path, "final_bins", "semibin_refined_bins.tsv.072_sub.fna")
            with open(source_genome2_path, "w") as f:
                f.write(">contigB\nTGCA\n")

            bin_info = pl.DataFrame([
                    [
                        "metabat2_refined_bins.tsv.005",
                        "p__Bacteroidetes (UID2605)", "350", "314", "208", "86.44", "1.8", "12.5", "1705643", "0", "315", "315", "6824", "6824", "5414", "5414", "25260", "25260", "59.9", "3.41", "90.22", "11",
                        "1669", "37", "269", "8", "0", "0", "0",
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "semibin_refined_bins.tsv.072_sub",
                        "o__Clostridiales (UID1212)", "172", "257", "149", "44.16", "1.03", "33.33", "980304", "0", "399", "399", "2442", "2442", "2456", "2456", "7735", "7735", "60.2", "2.97", "90.94", "11",
                        "1244", "143", "111", "3", "0", "0", "0",
                        "36.49", "0.24", "Neural Network (Specific Model)", "11", "0.909", "2442", "239.274115755627", "980304", "0.6", "1244", "399", "7735", ""
                    ],
                ], orient="row", schema=BIN_INFO_COLUMNS)
            (
                bin_info
                .write_csv(os.path.join(source_bins_path, "bin_info.tsv"), separator="\t")
            )

            cmd = (
                f"binchicken iterate "
                f"--coassemble-output prior_coassemble "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertIn("genomes", config)
            self.assertTrue(any("binchicken_coassembly_0.1.fna" in p for p in config["genomes"].values()))
            self.assertFalse(any("binchicken_coassembly_0.2.fna" in p for p in config["genomes"].values()))


if __name__ == '__main__':
    unittest.main()
