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
GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassembly_0")])

SAMPLE_READ_SIZE = os.path.join(MOCK_COASSEMBLE, "read_size2.csv")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "pipe", "sample_3_read.otu_table.tsv"),
    ])
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "summarise", "bins_summarised.otu_table2.tsv")

class Tests(unittest.TestCase):
    def test_iterate(self):
        with in_tempdir():
            cmd = (
                f"cockatoo iterate "
                f"--iteration 0 "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            new_bin_0_path = os.path.join("test", "recovered_bins", "iteration_0-coassembly_0-0.fna")
            self.assertTrue(os.path.exists(new_bin_0_path))
            with open(new_bin_0_path) as f:
                file = f.read()
                self.assertTrue(">k141_1363016" in file)
                self.assertTrue(">k141_177318" not in file)

            new_bin_1_path = os.path.join("test", "recovered_bins", "iteration_0-coassembly_0-1.fna")
            self.assertTrue(os.path.exists(new_bin_1_path))

            bin_provenance_path = os.path.join("test", "recovered_bins", "bin_provenance.tsv")
            self.assertTrue(os.path.exists(bin_provenance_path))

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

    def test_iterate_default_config(self):
        with in_tempdir():
            cmd = (
                f"cockatoo iterate "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--genome-transcripts {GENOME_TRANSCRIPTS} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["taxa_of_interest"], "")
            self.assertEqual(config["assemble_unmapped"], False)
            self.assertEqual(config["aviary_threads"], 16)
            self.assertEqual(config["aviary_memory"], 250)
            self.assertEqual(config["checkm_version"], 2)
            self.assertEqual(config["min_completeness"], 70)
            self.assertEqual(config["max_contamination"], 10)


if __name__ == '__main__':
    unittest.main()
