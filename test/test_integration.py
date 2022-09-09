#!/usr/bin/env python3

import unittest
import os
import sys
from bird_tool_utils import in_tempdir
import extern

path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','co.py')
path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

class Tests(unittest.TestCase):
    def test_coassemble(self):
        with in_tempdir():
            cmd = f"{path_to_script} coassemble --reads {path_to_data}/reads.fq --output test --threads 1"
            output = extern.run(cmd)

            expected = []
            self.assertEqual(expected, output)

    def test_coassemble_archive_otu_table(self):
        with in_tempdir():
            cmd = f"{path_to_script} coassemble --reads {path_to_data}/reads.fq --singlem-archive {path_to_data}/archive.json.gz --output test --threads 1"
            output = extern.run(cmd)

            expected = []
            self.assertEqual(expected, output)

    def test_evaluate(self):
        with in_tempdir():
            cmd = f"{path_to_script} evaluate --coassemble-output {path_to_data}/mock_coassemble --output test --threads 1"
            output = extern.run(cmd)

            expected = []
            self.assertEqual(expected, output)

if __name__ == '__main__':
    unittest.main()
