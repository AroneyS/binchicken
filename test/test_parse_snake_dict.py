#!/usr/bin/env python3

import unittest
from binchicken.common import parse_snake_dict

class Tests(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
