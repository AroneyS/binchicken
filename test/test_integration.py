#!/usr/bin/env python3

import unittest
import os
import sys
from bird_tool_utils import in_tempdir
import extern
from snakemake.io import load_configfile
from collections import OrderedDict

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

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
SAMPLE_SINGLEM = os.path.join(path_to_data, "sample_1.otu_table.tsv")
GENOMES = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
MOCK_CLUSTER = os.path.join(path_to_data, "mock_cluster")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_CLUSTER, "coassembly_0")])

class Tests(unittest.TestCase):
    def test_cluster(self):
        with in_tempdir():
            cmd = (
                f"cockatoo cluster "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--genome-transcripts {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "cluster", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "cluster", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "cluster", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_weight",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1,sample_2",
                        "2",
                        "2",
                        "2",
                        "0",
                        "sample_1,sample_2",
                        "coassembly_0"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_cluster_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo cluster "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--genome-transcripts {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["appraise_sequence_identity"], 0.89)
            self.assertEqual(config["max_coassembly_size"], 50)
            self.assertEqual(config["max_coassembly_samples"], 5)
            self.assertEqual(config["min_coassembly_coverage"], 10)
            self.assertEqual(config["max_recovery_samples"], 20)
            self.assertEqual(config["taxa_of_interest"], "")

    def test_coassemble(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--assemble-unmapped "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output $OUTPUT_DIR/coassembly_0/assemble -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_assemble.log ",
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
                        "aviary recover --assembly $OUTPUT_DIR/coassembly_0/assemble/assembly/final_contigs.fasta -1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output $OUTPUT_DIR/coassembly_0/recover -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_recover.log ",
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["assemble_unmapped"], False)

    def test_evaluate(self):
        with in_tempdir():
            cmd = (
                f"cockatoo evaluate "
                f"--cluster-output {MOCK_CLUSTER} "
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
                        "missed",
                        "recovered",
                        "total",
                        "recovered_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "1",
                        "1",
                        "2",
                        "50",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "0",
                        "1",
                        "1",
                        "100",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "1",
                        "1",
                        "2",
                        "50",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

    def test_evaluate_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo evaluate "
                f"--cluster-output {MOCK_CLUSTER} "
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

if __name__ == '__main__':
    unittest.main()
