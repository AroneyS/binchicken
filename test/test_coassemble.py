#!/usr/bin/env python3

import unittest
import os
import gzip
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

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2_protein.fna"),
    ])
MOCK_CLUSTER = os.path.join(path_to_data, "mock_cluster")

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

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

            test_dir = os.path.abspath("test")

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:118:HTHTVDSXX:1:1249:16740:14105" not in file)

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassembly_0", "assemble"),
                        "-n 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "logs", "coassembly_0_assemble.log"),
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
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassembly_0", "recover"),
                        "-n 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_abstract_options(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--assemble-unmapped "
                f"--genomes {GENOMES} "
                f"--abstract-options "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            test_dir = os.path.abspath("test")

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output $OUTPUT_DIR/coassembly_0/assemble -n $CPUS -m $MEMORY &> $OUTPUT_DIR/logs/coassembly_0_assemble.log "
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
                        "aviary recover --assembly $ASSEMBLE_OUTPUT/coassembly_0/assemble/assembly/final_contigs.fasta -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
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

            test_dir = os.path.abspath("test")

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
                        "--output", os.path.join(test_dir, "coassemble", "coassembly_0", "assemble"),
                        "-n 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "logs", "coassembly_0_assemble.log"),
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
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.1.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.1.fq"),
                        "-2",
                        os.path.join(os.path.abspath(path_to_data), "sample_1.2.fq"),
                        os.path.join(os.path.abspath(path_to_data), "sample_2.2.fq"),
                        "--output", os.path.join(test_dir, "coassemble", "coassembly_0", "recover"),
                        "-n 16 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_coassemble_genome_trim(self):
        with in_tempdir():
            extern.run(f"cp -R {MOCK_CLUSTER} mock_cluster")
            extern.run(f"cp {MOCK_CLUSTER}/appraise/binned.otu_table2.tsv mock_cluster/appraise/binned.otu_table.tsv")

            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output mock_cluster "
                f"--assemble-unmapped "
                f"--genomes {TWO_GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" finish_mapping\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            genomes_to_map_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertTrue(os.path.exists(genomes_to_map_path))
            with open(genomes_to_map_path) as f:
                lines = f.read()
                self.assertTrue("JABDGE010000038.1_13" in lines)
                self.assertTrue("TEST_CONTIG" not in lines)

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
            self.assertEqual(config["aviary_threads"], 16)
            self.assertEqual(config["aviary_memory"], 250)
            self.assertEqual(config["assemble_unmapped"], False)

    def test_coassemble_file_of_paths(self):
        with in_tempdir():
            write_string_to_file(SAMPLE_READS_FORWARD, "sample_reads_forward")
            write_string_to_file(SAMPLE_READS_REVERSE, "sample_reads_reverse")
            write_string_to_file(GENOMES, "genomes")

            cmd = (
                f"cockatoo coassemble "
                f"--forward-list sample_reads_forward "
                f"--reverse-list sample_reads_reverse "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--assemble-unmapped "
                f"--genomes-list genomes "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("collect_bins" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("coassemble_commands" in output)
            self.assertTrue("aviary_assemble" not in output)
            self.assertTrue("aviary_recover" not in output)
            self.assertTrue("collate_coassemblies" not in output)

    def test_coassemble_run_aviary(self):
        with in_tempdir():
            cmd = (
                f"cockatoo coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--cluster-output {MOCK_CLUSTER} "
                f"--assemble-unmapped "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--run-aviary "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("collect_bins" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("coassemble_commands" not in output)
            self.assertTrue("aviary_assemble" in output)
            self.assertTrue("aviary_recover" in output)
            self.assertTrue("collate_coassemblies" in output)


if __name__ == '__main__':
    unittest.main()
