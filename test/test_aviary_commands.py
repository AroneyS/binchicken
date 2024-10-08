#!/usr/bin/env python3

import unittest
import os
os.environ["POLARS_MAX_THREADS"] = "1"
import polars as pl
from polars.testing import assert_frame_equal
from binchicken.workflow.scripts.aviary_commands import pipeline

ELUSIVE_CLUSTERS_COLUMNS={
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
    }
COMMANDS_COLUMNS={
    "assemble": str,
    "recover": str,
}

reads_1 = {
    "1": "1_1.fq.gz",
    "2": "2_1.fq.gz",
    "3": "3_1.fq.gz",
    "4": "4_1.fq.gz",
    "5": "5_1.fq.gz",
    "6": "6_1.fq.gz",
}
reads_2 = {
    "1": "1_2.fq.gz",
    "2": "2_2.fq.gz",
    "3": "3_2.fq.gz",
    "4": "4_2.fq.gz",
    "5": "5_2.fq.gz",
    "6": "6_2.fq.gz",
}

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtypes=False)

    def test_aviary_commands(self):
        elusive_clusters = pl.DataFrame([
            ["1,2,3", 3, 2, 2000, "1,2,3,4,5,6", "coassembly_0"],
            ["4,5,6", 3, 1, 3000, "4,5,6", "coassembly_2"],
        ], orient="row", schema=ELUSIVE_CLUSTERS_COLUMNS)
        output_dir = "test_output_dir"
        assemble_threads = 10
        assemble_memory = 50
        recover_threads = 5
        recover_memory = 25


        expected_commands = pl.DataFrame([
            [
                f"aviary assemble --coassemble "
                f"-1 1_1.fq.gz 2_1.fq.gz 3_1.fq.gz "
                f"-2 1_2.fq.gz 2_2.fq.gz 3_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_0/assemble "
                f"-n {assemble_threads} -t {assemble_threads} -m {assemble_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_0_assemble.log ",

                f"aviary recover --assembly {output_dir}/coassemble/coassembly_0/assemble/assembly/final_contigs.fasta "
                f"-1 1_1.fq.gz 2_1.fq.gz 3_1.fq.gz 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 1_2.fq.gz 2_2.fq.gz 3_2.fq.gz 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_0/recover "
                f"-n {recover_threads} -t {recover_threads} -m {recover_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_0_recover.log ",
            ],
            [
                f"aviary assemble --coassemble "
                f"-1 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_2/assemble "
                f"-n {assemble_threads} -t {assemble_threads} -m {assemble_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_2_assemble.log ",

                f"aviary recover --assembly {output_dir}/coassemble/coassembly_2/assemble/assembly/final_contigs.fasta "
                f"-1 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_2/recover "
                f"-n {recover_threads} -t {recover_threads} -m {recover_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_2_recover.log ",
            ],
        ], orient="row", schema=COMMANDS_COLUMNS)

        observed_commands = pipeline(elusive_clusters, reads_1, reads_2, output_dir, assemble_threads, assemble_memory, recover_threads, recover_memory)
        self.assertDataFrameEqual(expected_commands, observed_commands)

    def test_aviary_commands_fast(self):
        elusive_clusters = pl.DataFrame([
            ["1,2,3", 3, 2, 2000, "1,2,3,4,5,6", "coassembly_0"],
            ["4,5,6", 3, 1, 3000, "4,5,6", "coassembly_2"],
        ], orient="row", schema=ELUSIVE_CLUSTERS_COLUMNS)
        output_dir = "test_output_dir"
        assemble_threads = 10
        assemble_memory = 50
        recover_threads = 5
        recover_memory = 25


        expected_commands = pl.DataFrame([
            [
                f"aviary assemble --coassemble "
                f"-1 1_1.fq.gz 2_1.fq.gz 3_1.fq.gz "
                f"-2 1_2.fq.gz 2_2.fq.gz 3_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_0/assemble "
                f"-n {assemble_threads} -t {assemble_threads} -m {assemble_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_0_assemble.log ",

                f"aviary recover --assembly {output_dir}/coassemble/coassembly_0/assemble/assembly/final_contigs.fasta "
                f"-1 1_1.fq.gz 2_1.fq.gz 3_1.fq.gz 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 1_2.fq.gz 2_2.fq.gz 3_2.fq.gz 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_0/recover "
                f"--binning-only "
                f"--refinery-max-iterations 0 "
                f"-n {recover_threads} -t {recover_threads} -m {recover_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_0_recover.log ",
            ],
            [
                f"aviary assemble --coassemble "
                f"-1 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_2/assemble "
                f"-n {assemble_threads} -t {assemble_threads} -m {assemble_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_2_assemble.log ",

                f"aviary recover --assembly {output_dir}/coassemble/coassembly_2/assemble/assembly/final_contigs.fasta "
                f"-1 4_1.fq.gz 5_1.fq.gz 6_1.fq.gz "
                f"-2 4_2.fq.gz 5_2.fq.gz 6_2.fq.gz "
                f"--output {output_dir}/coassemble/coassembly_2/recover "
                f"--binning-only "
                f"--refinery-max-iterations 0 "
                f"-n {recover_threads} -t {recover_threads} -m {recover_memory} --skip-qc "
                f"&> {output_dir}/coassemble/logs/coassembly_2_recover.log ",
            ],
        ], orient="row", schema=COMMANDS_COLUMNS)

        observed_commands = pipeline(elusive_clusters, reads_1, reads_2, output_dir, assemble_threads, assemble_memory, recover_threads, recover_memory, fast=True)
        self.assertDataFrameEqual(expected_commands, observed_commands)


if __name__ == '__main__':
    unittest.main()
