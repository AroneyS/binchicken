#!/usr/bin/env python3

import unittest
import os
import sys
from bird_tool_utils import in_tempdir
import extern

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

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
SAMPLE_SINGLEM = os.path.join(path_to_data, "sample_1.otu_table.tsv")
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")

class Tests(unittest.TestCase):
    def test_coassemble(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--output test"
            )
            output = extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))


            cluster_path = os.path.join("test", "coassembly", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))

    def test_coassemble_archive_otu_table(self):
        with in_tempdir():
            cmd = f"cockatoo coassemble --reads {SAMPLE_READS_FORWARD} --sample-singlem {SAMPLE_SINGLEM} --output test --dryrun"
            output = extern.run(cmd)

            expected = []
            self.assertEqual(expected, output)

    def test_evaluate(self):
        with in_tempdir():
            cmd = f"cockatoo evaluate --coassemble-output {MOCK_COASSEMBLE} --output test"
            output = extern.run(cmd)

            expected = []
            self.assertEqual(expected, output)

if __name__ == '__main__':
    unittest.main()
