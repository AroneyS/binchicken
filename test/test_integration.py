#!/usr/bin/env python3

import unittest
import os
import sys
from bird_tool_utils import in_tempdir
import extern

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
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_COASSEMBLIES = ' '.join([os.path.join(path_to_data, "mock_coassemble", "coassembly_0")])

class Tests(unittest.TestCase):
    def test_coassemble(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            edges_path = os.path.join("test", "coassembly", "target", "targets.tsv")
            self.assertTrue(os.path.exists(edges_path))

            edges_path = os.path.join("test", "coassembly", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(edges_path))

            cluster_path = os.path.join("test", "coassembly", "target", "elusive_clusters.tsv")
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

    def test_evaluate(self):
        with in_tempdir():
            cmd = (
                f"cockatoo evaluate "
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

if __name__ == '__main__':
    unittest.main()
