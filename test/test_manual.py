#!/usr/bin/env python3

import unittest
import os
import gzip
from bird_tool_utils import in_tempdir
import extern

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

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
APPRAISE_BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
APPRAISE_UNBINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "unbinned.otu_table.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters.tsv")
ELUSIVE_CLUSTERS_TWO = os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_two.tsv')


class Tests(unittest.TestCase):
    def test_update_sra_download_real(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--appraise-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
                f"--appraise-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
                f"--elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sra_1_path = os.path.join("test", "coassemble", "sra", "SRR8334323_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_1_path))
            with gzip.open(sra_1_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR8334323.1 HS2:487:H80UEADXX:1:1101:1148:1986/1" in file)
                self.assertTrue("@SRR8334323.2 HS2:487:H80UEADXX:1:1101:1148:1986/2" not in file)

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334323_2.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_1.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_2.fastq.gz")))

    def test_update_aviary_run_real(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--run-aviary "
                f"--cores 32 "
                f"--aviary-cores 32 "
                f"--aviary-memory 500 "
                f"--aviary-gtdbtk-dir /work/microbiome/db/gtdb/gtdb_release207_v2 "
                f"--aviary-checkm2-dir /work/microbiome/db/CheckM2_database "
                f"--genomes {GENOMES} "
                f"--appraise-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
                f"--appraise-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
                f"--elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" --config test=True\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            sra_1_path = os.path.join("test", "coassemble", "sra", "SRR8334323_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_1_path))
            with gzip.open(sra_1_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR8334323.1 1/1" in file)
                self.assertTrue("@SRR8334323.1 1/2" not in file)

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334323_2.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_1.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_2.fastq.gz")))

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fna")))

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "coassemble", "coassembly_0", "recover", "bins", "checkm_minimal.tsv")))

    def test_update_specific_coassembly_sra(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--appraise-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
                f"--appraise-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
                f"--elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra_two.tsv')} "
                f"--coassemblies coassembly_0 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            output = extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334323_1.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334323_2.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_1.fastq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334324_2.fastq.gz")))

            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334320_1.fastq.gz")))
            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334320_2.fastq.gz")))
            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334321_1.fastq.gz")))
            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "sra", "SRR8334321_2.fastq.gz")))


if __name__ == '__main__':
    unittest.main()
