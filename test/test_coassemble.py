#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
from snakemake.io import load_configfile

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

GENOMES = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
MOCK_CLUSTER = os.path.join(path_to_data, "mock_cluster")

class Tests(unittest.TestCase):
    def test_coassemble(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--assemble-unmapped "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")),
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.1.fq.gz")),
                        "-2",
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.2.fq.gz")),
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.2.fq.gz")),
                        "--output $OUTPUT_DIR/coassembly_0/assemble -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_assemble.log ",
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
                        "aviary recover --assembly $OUTPUT_DIR/coassembly_0/assemble/assembly/final_contigs.fasta -1",
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")),
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.1.fq.gz")),
                        "-2",
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.2.fq.gz")),
                        os.path.abspath(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.2.fq.gz")),
                        "--output $OUTPUT_DIR/coassembly_0/recover -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_recover.log ",
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_no_mapping(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output $OUTPUT_DIR/coassembly_0/assemble -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_assemble.log ",
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
                        "aviary recover --assembly $OUTPUT_DIR/coassembly_0/assemble/assembly/final_contigs.fasta -1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output $OUTPUT_DIR/coassembly_0/recover -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_recover.log ",
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["assemble_unmapped"], False)


if __name__ == '__main__':
    unittest.main()
