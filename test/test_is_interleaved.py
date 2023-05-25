#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from ibis.workflow.scripts.is_interleaved import pipeline, is_interleaved

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','ibis','workflow','scripts','is_interleaved.py')

path_to_sra = os.path.join(path_to_data, "sra")
SRA_INTERLEAVED = os.path.join(path_to_sra, "SRR3309137.fastq")
SRA_INTERLEAVED_GZ = os.path.join(path_to_sra, "SRR3309137.fastq.gz")
SRA_ODD = os.path.join(path_to_sra, "SRR3309137_odd.fastq")
SRA_ODD_GZ = os.path.join(path_to_sra, "SRR3309137_odd.fastq.gz")
SRA_MISMATCHED = os.path.join(path_to_sra, "SRR3309137_mismatched.fastq")
SRA_MISMATCHED_GZ = os.path.join(path_to_sra, "SRR3309137_mismatched.fastq.gz")
SRA_EMPTY = os.path.join(path_to_sra, "SRR3309137_empty.fastq")

READS_COLUMNS={
    "read": str,
    "number": int,
    }

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False)

    def test_is_interleaved(self):
        reads = pl.DataFrame([
            ["@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 1],
            ["@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 2],
            ["@SRR3309137.3 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 3],
            ["@SRR3309137.4 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 4],
            ["@SRR3309137.5 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 5],
            ["@SRR3309137.6 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 6],
            ["@SRR3309137.7 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 7],
            ["@SRR3309137.8 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 8],
            ["@SRR3309137.9 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 9],
            ["@SRR3309137.10 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 10],
            ["@SRR3309137.104315991 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315991],
            ["@SRR3309137.104315992 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315992],
            ["@SRR3309137.104315993 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315993],
            ["@SRR3309137.104315994 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315994],
            ["@SRR3309137.104315995 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315995],
            ["@SRR3309137.104315996 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315996],
            ["@SRR3309137.104315997 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315997],
            ["@SRR3309137.104315998 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315998],
            ["@SRR3309137.104315999 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104315999],
            ["@SRR3309137.104316000 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104316000],
        ], schema=READS_COLUMNS)
        total_count = 104316000

        observed_outcome, observed_reason = is_interleaved(reads, total_count)
        self.assertEqual("Duplicate read names in consecutive reads with even readcount", observed_reason)
        self.assertEqual(True, observed_outcome)

    def test_is_interleaved_odd_count(self):
        reads = pl.DataFrame([
            ["@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 1],
            ["@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 2],
            ["@SRR3309137.3 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 3],
            ["@SRR3309137.4 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 4],
            ["@SRR3309137.5 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 5],
            ["@SRR3309137.6 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 6],
            ["@SRR3309137.7 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 7],
            ["@SRR3309137.8 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 8],
            ["@SRR3309137.9 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 9],
            ["@SRR3309137.10 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 10],
            ["@SRR3309137.104315991 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315991],
            ["@SRR3309137.104315992 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315992],
            ["@SRR3309137.104315993 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315993],
            ["@SRR3309137.104315994 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315994],
            ["@SRR3309137.104315995 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315995],
            ["@SRR3309137.104315996 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315996],
            ["@SRR3309137.104315997 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315997],
            ["@SRR3309137.104315998 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315998],
            ["@SRR3309137.104315999 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104315999],
            ["@SRR3309137.104316000 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104316000],
        ], schema=READS_COLUMNS)
        total_count = 104316001

        observed_outcome, observed_reason = is_interleaved(reads, total_count)
        self.assertEqual("Odd readcount", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_mismatched(self):
        reads = pl.DataFrame([
            ["@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 1],
            ["@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:9999:9999/1", 2],
            ["@SRR3309137.3 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 3],
            ["@SRR3309137.4 HISEQ06:195:D1DRHACXX:5:1101:9999:9999/1", 4],
            ["@SRR3309137.5 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 5],
            ["@SRR3309137.6 HISEQ06:195:D1DRHACXX:5:1101:9999:9999/1", 6],
            ["@SRR3309137.7 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 7],
            ["@SRR3309137.8 HISEQ06:195:D1DRHACXX:5:1101:9999:9999/1", 8],
            ["@SRR3309137.9 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 9],
            ["@SRR3309137.10 HISEQ06:195:D1DRHACXX:5:1101:9999:9999/1", 10],
            ["@SRR3309137.104315991 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315991],
            ["@SRR3309137.104315992 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104315992],
            ["@SRR3309137.104315993 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315993],
            ["@SRR3309137.104315994 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104315994],
            ["@SRR3309137.104315995 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315995],
            ["@SRR3309137.104315996 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104315996],
            ["@SRR3309137.104315997 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315997],
            ["@SRR3309137.104315998 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104315998],
            ["@SRR3309137.104315999 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104315999],
            ["@SRR3309137.104316000 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104316000],
        ], schema=READS_COLUMNS)
        total_count = 104316000

        observed_outcome, observed_reason = is_interleaved(reads, total_count)
        self.assertEqual("Consecutive reads do not match (10/10)", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_mismatched_some(self):
        reads = pl.DataFrame([
            ["@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 1],
            ["@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1", 2],
            ["@SRR3309137.3 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 3],
            ["@SRR3309137.4 HISEQ06:195:D1DRHACXX:5:1101:1544:2238/1", 4],
            ["@SRR3309137.5 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 5],
            ["@SRR3309137.6 HISEQ06:195:D1DRHACXX:5:1101:1834:2245/1", 6],
            ["@SRR3309137.7 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 7],
            ["@SRR3309137.8 HISEQ06:195:D1DRHACXX:5:1101:2382:2240/1", 8],
            ["@SRR3309137.9 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 9],
            ["@SRR3309137.10 HISEQ06:195:D1DRHACXX:5:1101:2690:2245/1", 10],
            ["@SRR3309137.104315991 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315991],
            ["@SRR3309137.104315992 HISEQ06:195:D1DRHACXX:5:2308:21349:200557/1", 104315992],
            ["@SRR3309137.104315993 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315993],
            ["@SRR3309137.104315994 HISEQ06:195:D1DRHACXX:5:2308:21316:200559/1", 104315994],
            ["@SRR3309137.104315995 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315995],
            ["@SRR3309137.104315996 HISEQ06:195:D1DRHACXX:5:2308:21389:200636/1", 104315996],
            ["@SRR3309137.104315997 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315997],
            ["@SRR3309137.104315998 HISEQ06:195:D1DRHACXX:5:2308:21283:200652/1", 104315998],
            ["@SRR3309137.104315999 HISEQ06:195:D1DRHACXX:5:2308:15935:200751/1", 104315999],
            ["@SRR3309137.104316000 HISEQ06:195:D1DRHACXX:5:2308:99999:999999/1", 104316000],
        ], schema=READS_COLUMNS)
        total_count = 104316000

        observed_outcome, observed_reason = is_interleaved(reads, total_count)
        self.assertEqual("Consecutive reads do not match (1/10)", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline(self):
        observed_outcome, observed_reason = pipeline(SRA_INTERLEAVED, 5, 5)
        self.assertEqual("Duplicate read names in consecutive reads with even readcount", observed_reason)
        self.assertEqual(True, observed_outcome)

    def test_is_interleaved_pipeline_gz(self):
        observed_outcome, observed_reason = pipeline(SRA_INTERLEAVED_GZ, 5, 5)
        self.assertEqual("Duplicate read names in consecutive reads with even readcount", observed_reason)
        self.assertEqual(True, observed_outcome)

    def test_is_interleaved_pipeline_odd(self):
        observed_outcome, observed_reason = pipeline(SRA_ODD, 5, 5)
        self.assertEqual("Odd readcount", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline_odd_gz(self):
        observed_outcome, observed_reason = pipeline(SRA_ODD_GZ, 5, 5)
        self.assertEqual("Odd readcount", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline_mismatched(self):
        observed_outcome, observed_reason = pipeline(SRA_MISMATCHED, 5, 5)
        self.assertEqual("Consecutive reads do not match (1/10)", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline_mismatched_gz(self):
        observed_outcome, observed_reason = pipeline(SRA_MISMATCHED_GZ, 5, 5)
        self.assertEqual("Consecutive reads do not match (1/10)", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline_mismatched_shorter(self):
        observed_outcome, observed_reason = pipeline(SRA_MISMATCHED, 2, 3)
        self.assertEqual("Consecutive reads do not match (1/5)", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_pipeline_empty(self):
        observed_outcome, observed_reason = pipeline(SRA_EMPTY, 2, 3)
        self.assertEqual("Empty file", observed_reason)
        self.assertEqual(False, observed_outcome)

    def test_is_interleaved_cli(self):
        with in_tempdir():
            cmd = (
                f"python {path_to_script} "
                f"--input {SRA_INTERLEAVED} "
                f"--start-check-pairs 5 "
                f"--end-check-pairs 5 "
            )
            output = extern.run(cmd)

            output_text = output.strip().split("\t")
            self.assertEqual("True", output_text[0])
            self.assertEqual("Duplicate read names in consecutive reads with even readcount", output_text[1])


if __name__ == '__main__':
    unittest.main()
