#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from binchicken.workflow.scripts.summarise_coassemblies import processing

ELUSIVE_CLUSTERS_COLUMNS={
    "coassembly": str,
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    }

READ_SIZE_COLUMNS={
    "sample": str,
    "read_size": int,
    }

SUMMARY_COLUMNS={
    "coassembly": str,
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "unmapped_size": int,
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False, check_row_order=False)

    def test_summarise_coassemblies(self):
        elusive_clusters = pl.DataFrame([
            ["coassembly_1", "sample_1,sample_2", 2, 200, 1000],
            ["coassembly_2", "sample_1,sample_3", 2, 100, 2000],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        read_size = pl.DataFrame([
            ["sample_1", 100],
            ["sample_2", 200],
            ["sample_3", 300],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([
            ["coassembly_1", "sample_1,sample_2", 2, 200, 1000, 300],
            ["coassembly_2", "sample_1,sample_3", 2, 100, 2000, 400],
        ], schema=SUMMARY_COLUMNS)

        observed = processing(elusive_clusters, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_summarise_coassemblies_empty_input(self):
        elusive_clusters = pl.DataFrame([], schema=ELUSIVE_CLUSTERS_COLUMNS)
        read_size = pl.DataFrame([], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([], schema=SUMMARY_COLUMNS)

        observed = processing(elusive_clusters, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_summarise_coassemblies_empty_elusive_clusters(self):
        elusive_clusters = pl.DataFrame([], schema=ELUSIVE_CLUSTERS_COLUMNS)
        read_size = pl.DataFrame([
            ["sample_1", 100],
            ["sample_2", 200],
            ["sample_3", 300],
        ], schema=READ_SIZE_COLUMNS)

        expected = pl.DataFrame([], schema=SUMMARY_COLUMNS)

        observed = processing(elusive_clusters, read_size)
        self.assertDataFrameEqual(expected, observed)

    def test_summarise_coassemblies_empty_read_size(self):
        elusive_clusters = pl.DataFrame([
            ["coassembly_1", "sample_1,sample_2", 2, 200, 1000],
            ["coassembly_2", "sample_1,sample_3", 2, 100, 2000],
        ], schema=ELUSIVE_CLUSTERS_COLUMNS)

        read_size = None

        expected = pl.DataFrame([
            ["coassembly_1", "sample_1,sample_2", 2, 200, 1000],
            ["coassembly_2", "sample_1,sample_3", 2, 100, 2000],
        ], schema={k:v for k,v in SUMMARY_COLUMNS.items() if k != "unmapped_size"})

        observed = processing(elusive_clusters, read_size)
        self.assertDataFrameEqual(expected, observed)


if __name__ == '__main__':
    unittest.main()
