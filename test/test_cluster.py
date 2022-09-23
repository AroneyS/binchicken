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
GENOMES = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
MOCK_CLUSTER = os.path.join(path_to_data, "mock_cluster")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_CLUSTER, "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_CLUSTER, "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_CLUSTER, "pipe", "sample_3_read.otu_table.tsv"),
    ])
SAMPLE_READ_SIZE = os.path.join(MOCK_CLUSTER, "read_size2.csv")
GENOME_SINGLEM = os.path.join(MOCK_CLUSTER, "summarise", "bins_summarised.otu_table2.tsv")
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

            read_size_path = os.path.join("test", "cluster", "read_size.csv")
            self.assertTrue(os.path.exists(read_size_path))
            expected = "\n".join(
                [
                    ",".join(["sample_1", "1661"]),
                    ",".join(["sample_2", "1208"]),
                    ",".join(["sample_3", "1208"]),
                    ""
                ]
            )
            with open(read_size_path) as f:
                self.assertEqual(expected, f.read())

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
                        "2869",
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
            self.assertEqual(config["taxa_of_interest"], "")

    def test_cluster_singlem_inputs(self):
        with in_tempdir():
            cmd = (
                f"cockatoo cluster "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("singlem_pipe_bins" not in output)
            self.assertTrue("singlem_summarise_bins" not in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("target_elusive" in output)
            self.assertTrue("cluster_graph" in output)

    def test_cluster_taxa_of_interest(self):
        with in_tempdir():
            cmd = (
                f"cockatoo cluster "
                f"--sample-singlem {SAMPLE_SINGLEM} "
                f"--sample-read-size {SAMPLE_READ_SIZE} "
                f"--genome-singlem {GENOME_SINGLEM} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--taxa-of-interest \"p__Actinobacteriota\" "
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
                        "sample_1,sample_3",
                        "2",
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


if __name__ == '__main__':
    unittest.main()
