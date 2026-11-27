#!/usr/bin/env python3

import unittest
from binchicken.common import split_bin_name, fasta_to_name, parse_snake_dict, bin_coassembly

class Tests(unittest.TestCase):
    def test_fasta_to_name(self):
        filename1 = "data/sample.1.fq.gz"

        expected = "sample.1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_R(self):
        filename1 = "data/sample_R1.fq.gz"

        expected = "sample_R1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_1(self):
        filename1 = "data/sample_1.fq.gz"

        expected = "sample_1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_fna_gz(self):
        filename1 = "data/sample_1.fna.gz"

        expected = "sample_1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_fna(self):
        filename1 = "data/sample_1.fna"

        expected = "sample_1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_fa(self):
        filename1 = "data/sample_1.fa"

        expected = "sample_1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_fasta_to_name_fasta(self):
        filename1 = "data/sample_1.fasta"

        expected = "sample_1"
        observed = fasta_to_name(filename1)
        self.assertEqual(expected, observed)

    def test_parse_snake_dict(self):
        snake_dict = '{sample_1: /hpccs01/home/aroneys/src/binchicken/test/data/sample_1.1.fq, sample_2: /hpccs01/home/aroneys/src/binchicken/test/data/sample_2.1.fq, sample_3: /hpccs01/home/aroneys/src/binchicken/test/data/sample_3.1.fq, sample_5: /hpccs01/home/aroneys/src/binchicken/test/data/sample_5.1.fq}'

        expected = {
            "sample_1": "/hpccs01/home/aroneys/src/binchicken/test/data/sample_1.1.fq",
            "sample_2": "/hpccs01/home/aroneys/src/binchicken/test/data/sample_2.1.fq",
            "sample_3": "/hpccs01/home/aroneys/src/binchicken/test/data/sample_3.1.fq",
            "sample_5": "/hpccs01/home/aroneys/src/binchicken/test/data/sample_5.1.fq",
            }
        observed = parse_snake_dict(snake_dict)
        self.assertEqual(expected, observed)

    def test_parse_snake_dict_empty(self):
        snake_dict = '{}'

        expected = {}
        observed = parse_snake_dict(snake_dict)
        self.assertEqual(expected, observed)

    def test_parse_snake_dict_genome(self):
        snake_dict = '{GB_GCA_013286235.1: /home/runner/work/binchicken/binchicken/test/data/GB_GCA_013286235.1.fna}'

        expected = {
            "GB_GCA_013286235.1": "/home/runner/work/binchicken/binchicken/test/data/GB_GCA_013286235.1.fna",
            }
        observed = parse_snake_dict(snake_dict)
        self.assertEqual(expected, observed)

    def test_parse_snake_dict_special(self):
        snake_dict = '{sample-1: /path/to/file-1.fq, sample_2: /path/to/file_2.fq}'

        expected = {
            "sample-1": "/path/to/file-1.fq",
            "sample_2": "/path/to/file_2.fq"
        }
        observed = parse_snake_dict(snake_dict)
        self.assertEqual(expected, observed)

    def test_parse_snake_dict_numeric(self):
        snake_dict = '{1: /path/to/file1.fq, 2: /path/to/file2.fq}'

        expected = {
            "1": "/path/to/file1.fq",
            "2": "/path/to/file2.fq"
        }
        observed = parse_snake_dict(snake_dict)
        self.assertEqual(expected, observed)


    def test_bin_coassembly(self):
        cases = {
            "binchicken_coassembly_1.3": "coassembly_1",
            "binchicken_coassembly_10.12": "coassembly_10",
            "binchicken_it2_coassembly_0.1": "coassembly_0",
            "binchicken_sample_1.1": "sample_1",
            "binchicken_SRR123.4.2": "SRR123.4",
            "coassembly_1-0": "coassembly_1",
            "coassembly_10-3": "coassembly_10",
            "sample_1-0": "sample_1",
            "coassembly_0-GB_GCA_013286235.1": "coassembly_0",
        }
        for name, expected in cases.items():
            self.assertEqual(expected, bin_coassembly(name), name)


    def test_bin_coassembly_known_names(self):
        # With known coassembly names, names are resolved against them rather than parsed blindly
        known = ["it2_sample", "sample-a", "sample-a-b", "sample-b", "coassembly_0"]
        cases = {
            "binchicken_it2_sample.1": "it2_sample",
            "binchicken_it1_it2_sample.1": "it2_sample",
            "binchicken_it1_coassembly_0.3": "coassembly_0",
            "sample-a-0": "sample-a",
            "sample-b-0": "sample-b",
            "sample-a-b-0": "sample-a-b",
            "coassembly_0-genome-1": "coassembly_0",
        }
        for name, expected in cases.items():
            self.assertEqual(expected, bin_coassembly(name, known), name)

    def test_split_bin_name(self):
        known = ["sample-a", "it2_sample"]
        cases = {
            # Aviary/new-genome keys keep the remaining genome suffix
            "sample-a-0_transcripts": ("sample-a", "0_transcripts"),
            "coassembly_0-genome-1": ("coassembly_0", "genome-1"),
            # Recovered names keep the whole name as the genome
            "binchicken_it2_sample.1_transcripts": ("it2_sample", "binchicken_it2_sample.1_transcripts"),
            "binchicken_it2_sample.1": ("it2_sample", "binchicken_it2_sample.1"),
        }
        for name, expected in cases.items():
            self.assertEqual(expected, split_bin_name(name, known), name)

if __name__ == '__main__':
    unittest.main()
