#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
import subprocess
import polars as pl
from polars.testing import assert_frame_equal

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
SAMPLE_READS_FORWARD_PRE = " ".join([SAMPLE_READS_FORWARD, os.path.join(path_to_data, "sample_5.1.fq")])
SAMPLE_READS_REVERSE_PRE = " ".join([SAMPLE_READS_REVERSE, os.path.join(path_to_data, "sample_5.2.fq")])

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
SAMPLE_SINGLEM_PRE = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_3_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe_precluster", "sample_5_read.otu_table.tsv"),
    ])

PRECLUSTER_DISTANCES = os.path.join(path_to_data, "sketch", "sample_distances.csv")

PRIOR_MISSING = os.path.join(path_to_data, "prior_missing.tsv")
PRIOR_EXTRA = os.path.join(path_to_data, "prior_extra.tsv")

class Tests(unittest.TestCase):
    def test_single(self):
        with in_tempdir():
            cmd = (
                f"binchicken single "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--coassembly-samples sample_1 sample_2 "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "4",
                        "4832",
                        "sample_1,sample_2,sample_3",
                        "sample_1"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "3",
                        "3926",
                        "sample_1,sample_2,sample_3",
                        "sample_2"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            test_dir = os.path.abspath("test")
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble --coassemble -1",
                        SAMPLE_READS_FORWARD.split(" ")[0],
                        "-2",
                        SAMPLE_READS_REVERSE.split(" ")[0],
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "sample_1", "assemble"),
                        "-n 64 -t 64 -m 500 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "sample_1_assemble.log"),
                        ""
                    ]),
                    " ".join([
                        "aviary assemble --coassemble -1",
                        SAMPLE_READS_FORWARD.split(" ")[1],
                        "-2",
                        SAMPLE_READS_REVERSE.split(" ")[1],
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "sample_2", "assemble"),
                        "-n 64 -t 64 -m 500 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "sample_2_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "sample_1", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        SAMPLE_READS_FORWARD,
                        "-2",
                        SAMPLE_READS_REVERSE,
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "sample_1", "recover"),
                        "--binning-only --refinery-max-iterations 0 "
                        "-n 32 -t 32 -m 250 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "sample_1_recover.log"),
                        ""
                    ]),
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "sample_2", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        SAMPLE_READS_FORWARD,
                        "-2",
                        SAMPLE_READS_REVERSE,
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "sample_2", "recover"),
                        "--binning-only --refinery-max-iterations 0 "
                        "-n 32 -t 32 -m 250 --skip-qc &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "sample_2_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

            summary_path = os.path.join("test", "coassemble", "summary.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                    ]),
                    "\t".join([
                        "sample_1",
                        "sample_1",
                        "1",
                        "4",
                        "4832",
                    ]),
                    "\t".join([
                        "sample_2",
                        "sample_2",
                        "1",
                        "3",
                        "3926",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

    def test_single_preclustered(self):
        with in_tempdir():
            cmd = (
                f"binchicken single "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--max-recovery-samples 2 "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_2,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_3,sample_5", "1779956245962588283,5034568815038442683"],
                ],
                schema = ["style", "cluster_size", "samples", "target_ids"],
                orient="row",
            )
            observed = pl.read_csv(elusive_edges_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False)

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_5",
                        "1",
                        "3",
                        "3624",
                        "sample_3,sample_5",
                        "sample_5"
                    ]),
                    "\t".join([
                        "sample_3",
                        "1",
                        "2",
                        "3624",
                        "sample_3,sample_5",
                        "sample_3"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "2",
                        "3926",
                        "sample_1,sample_2",
                        "sample_2"
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "2",
                        "4832",
                        "sample_1,sample_2",
                        "sample_1"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_single_preclustered_distances_provided(self):
        with in_tempdir():
            cmd = (
                f"binchicken single "
                f"--forward {SAMPLE_READS_FORWARD_PRE} "
                f"--reverse {SAMPLE_READS_REVERSE_PRE} "
                f"--sample-singlem {SAMPLE_SINGLEM_PRE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--precluster-distances {PRECLUSTER_DISTANCES} "
                f"--kmer-precluster always "
                f"--precluster-size 3 "
                f"--max-recovery-samples 2 "
                f"--output test "
            )
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')

            self.assertTrue("sketch_samples" not in output)
            self.assertTrue("distance_samples" not in output)
            self.assertTrue("provided_distances" in output)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sketch_path = os.path.join("test", "coassemble", "sketch", "samples.sig")
            self.assertFalse(os.path.exists(sketch_path))

            distances_path = os.path.join("test", "coassemble", "sketch", "samples.csv")
            self.assertTrue(os.path.exists(distances_path))

            elusive_edges_path = os.path.join("test", "coassemble", "target", "elusive_edges.tsv")
            self.assertTrue(os.path.exists(elusive_edges_path))
            expected = pl.DataFrame([
                    ["match", 2, "sample_1,sample_2", "2416810233479675551,6360827971060584481"],
                    ["match", 2, "sample_1,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_2,sample_5", "2416810233479675551"],
                    ["match", 2, "sample_3,sample_5", "1779956245962588283,5034568815038442683"],
                ],
                schema = ["style", "cluster_size", "samples", "target_ids"],
                orient="row",
            )
            observed = pl.read_csv(elusive_edges_path, separator="\t")
            assert_frame_equal(expected, observed, check_dtypes=False, check_row_order=False)

            cluster_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "samples",
                        "length",
                        "total_targets",
                        "total_size",
                        "recover_samples",
                        "coassembly",
                    ]),
                    "\t".join([
                        "sample_5",
                        "1",
                        "3",
                        "3624",
                        "sample_3,sample_5",
                        "sample_5"
                    ]),
                    "\t".join([
                        "sample_3",
                        "1",
                        "2",
                        "3624",
                        "sample_3,sample_5",
                        "sample_3"
                    ]),
                    "\t".join([
                        "sample_2",
                        "1",
                        "2",
                        "3926",
                        "sample_1,sample_2",
                        "sample_2"
                    ]),
                    "\t".join([
                        "sample_1",
                        "1",
                        "2",
                        "4832",
                        "sample_1,sample_2",
                        "sample_1"
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_single_assembly_missing(self):
        with in_tempdir():
            cmd = (
                f"binchicken single "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--coassembly-samples sample_1 sample_2 "
                f"--prior-assemblies {PRIOR_MISSING} "
                f"--output test "
            )

            with self.assertRaises(Exception) as context:
                _ = extern.run(cmd)

            self.assertTrue("Samples/coassemblies missing assemblies in prior assemblies: sample_2" in str(context.exception))

    def test_single_assembly_extra(self):
        with in_tempdir():
            cmd = (
                f"binchicken single "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--coassembly-samples sample_1 sample_2 "
                f"--prior-assemblies {PRIOR_EXTRA} "
                f"--output test "
            )

            with self.assertRaises(Exception) as context:
                _ = extern.run(cmd)

            self.assertTrue("Extra assemblies not matching any samples/coassemblies in prior assemblies: sample_3" in str(context.exception))

if __name__ == '__main__':
    unittest.main()
