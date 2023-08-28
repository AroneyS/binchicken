#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
import subprocess
from snakemake.io import load_configfile
import re

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
GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_UNBINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "unbinned.otu_table.tsv")
MOCK_BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0")])
MOCK_GENOMES = " ".join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_1.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_2.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_3.fna"),
])
MOCK_GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "summarise", "bins_summarised_mock.otu_table.tsv")
ELUSIVE_CLUSTERS = ' '.join([os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters_iterate.tsv")])

SAMPLE_READ_SIZE = os.path.join(MOCK_COASSEMBLE, "read_size2.csv")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_3_read.otu_table.tsv"),
    ])
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "summarise", "bins_summarised.otu_table2.tsv")

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

class Tests(unittest.TestCase):
    def test_iterate(self):
        with in_tempdir():
            cmd = (
                f"ibis iterate "
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
                f"--conda-prefix {path_to_conda} "
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

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
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

    def test_iterate_genome_input(self):
        with in_tempdir():
            cmd = (
                f"ibis iterate "
                f"--new-genomes {MOCK_GENOMES} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--prodigal-meta "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            new_bin_0_path = os.path.join("test", "recovered_bins", "bin_1.fna")
            self.assertTrue(os.path.exists(new_bin_0_path))
            with open(new_bin_0_path) as f:
                file = f.read()
                self.assertTrue(">k141_1363016" in file)
                self.assertTrue(">k141_177318" not in file)

            new_bin_1_path = os.path.join("test", "recovered_bins", "bin_2.fna")
            self.assertTrue(os.path.exists(new_bin_1_path))

            new_bin_2_path = os.path.join("test", "recovered_bins", "bin_3.fna")
            self.assertTrue(os.path.exists(new_bin_2_path))

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

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
                f"ibis iterate "
                f"--new-genomes-list genomes "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("count_bp_reads" in output)
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
                os.path.splitext(os.path.basename(g))[0]: g.replace(MOCK_COASSEMBLE + "/coassemble/coassembly_0/recover/bins/final_bins/", os.getcwd() + "/test/recovered_bins/")
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
                f"ibis iterate "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
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
            self.assertEqual(config["aviary_threads"], 64)
            self.assertEqual(config["aviary_memory"], 500)

            self.assertTrue("Evaluating bins using CheckM2 with completeness >= 70 and contamination <= 10" in output)

    def test_iterate_genome_singlem(self):
        with in_tempdir():
            cmd = (
                f"ibis iterate "
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
                f"--conda-prefix {path_to_conda} "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
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
                f"ibis iterate "
                f"--iteration 0 "
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
                f"--conda-prefix {path_to_conda} "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
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
                f"ibis iterate "
                f"--iteration 0 "
                f"--new-genomes {MOCK_GENOMES} "
                f"--new-genome-singlem {MOCK_GENOME_SINGLEM} "
                f"--elusive-clusters {ELUSIVE_CLUSTERS} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
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


if __name__ == '__main__':
    unittest.main()
