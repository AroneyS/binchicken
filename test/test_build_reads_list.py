#!/usr/bin/env python3

import unittest
import os

from binchicken.binchicken import build_reads_list

class Tests(unittest.TestCase):
    def test_build_reads_list(self):
        # Unsorted inputs that would mismatch without sorting
        forward = [
            os.path.join("reads", "XS-21-1_1.fastq"),
            os.path.join("reads", "XS-21-12_1.fastq"),
        ]
        reverse = [
            os.path.join("reads", "XS-21-12_2.fastq"),
            os.path.join("reads", "XS-21-1_2.fastq"),
        ]

        f_map, r_map = build_reads_list(forward, reverse, no_sample_sort=False)

        f_keys = list(f_map.keys())
        r_keys = list(r_map.keys())
        self.assertEqual(f_keys, r_keys)

    def test_build_reads_list_no_sample_sort(self):
        # Same inputs, but no_sample_sort preserves original zip order
        forward = [
            os.path.join("reads", "XS-21-1_1.fastq"),
            os.path.join("reads", "XS-21-12_1.fastq"),
        ]
        reverse = [
            os.path.join("reads", "XS-21-12_2.fastq"),
            os.path.join("reads", "XS-21-1_2.fastq"),
        ]

        # Expect duplicate basename error due to mismatched zip order
        with self.assertRaises(Exception) as ctx:
            build_reads_list(forward, reverse, no_sample_sort=True)
        self.assertIn("Duplicate basename: XS-21-1", str(ctx.exception))

if __name__ == "__main__":
    unittest.main()
