#!/usr/bin/env python3

import unittest
import os
import shutil
import gzip
import subprocess
from ruamel.yaml import YAML
import extern
from bird_tool_utils import in_tempdir
import pytest

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

SINGLEM_METAPACKAGE = "/work/microbiome/db/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"
GTDBTK_DB = "/work/microbiome/db/gtdb/gtdb_release207_v2"
CHECKM2_DB = "/work/microbiome/db/CheckM2_database"
METABULI_DB = "/work/microbiome/abisko/aroneys/db/metabuli"

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
SAMPLE_READS_FORWARD_EMPTY = " ".join([SAMPLE_READS_FORWARD, os.path.join(path_to_data, "sample_4.1.fq")])
SAMPLE_READS_REVERSE_EMPTY = " ".join([SAMPLE_READS_REVERSE, os.path.join(path_to_data, "sample_4.2.fq")])

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
SAMPLE_READ_SIZE = os.path.join(MOCK_COASSEMBLE, "coassemble", "read_size2.csv")
SAMPLE_SINGLEM = ' '.join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_1_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_2_read.otu_table.tsv"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "pipe", "sample_3_read.otu_table.tsv"),
    ])

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
GENOME_TRANSCRIPTS = ' '.join([os.path.join(path_to_data, "GB_GCA_013286235.1_protein.fna")])
GENOME_SINGLEM = os.path.join(MOCK_COASSEMBLE, "coassemble", "summarise", "bins_summarised.otu_table2.tsv")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble", "coassemble")
APPRAISE_BINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "binned.otu_table.tsv")
APPRAISE_UNBINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned.otu_table.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_clusters.tsv")
ELUSIVE_CLUSTERS_TWO = os.path.join(MOCK_COASSEMBLE, "coassemble", 'target', 'elusive_clusters_two.tsv')

MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
PRIOR_ASSEMBLY = os.path.join(path_to_data, "prior_assembly.tsv")
PRIOR_COASSEMBLY = os.path.join(path_to_data, "prior_coassembly.tsv")

@pytest.mark.expensive
class Tests(unittest.TestCase):
    def setup_output_dir(self, output_dir):
        try:
            shutil.rmtree(output_dir)
        except FileNotFoundError:
            pass
        os.makedirs(output_dir)

    def test_coassemble_sra_download_real(self):
        output_dir = os.path.join("example", "test_coassemble_sra_download_real")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken coassemble "
            f"--forward SRR8334323 SRR8334324 SRR6797127 SRR6797128 "
            f"--singlem-metapackage {SINGLEM_METAPACKAGE} "
            f"--sra "
            f"--cores 32 "
            f"--output {output_dir} "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR6797127_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR6797127_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR6797128_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR6797128_2.fastq.gz")))

        cluster_path = os.path.join(output_dir, "coassemble", "target", "elusive_clusters.tsv")
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
                    "SRR8334323,SRR8334324",
                    "2",
                    "29",
                    "5046179100",
                    "SRR8334323,SRR8334324",
                    "coassembly_0"
                ]),
                "\t".join([
                    "SRR6797127,SRR6797128",
                    "2",
                    "10",
                    "7675876600",
                    "SRR6797127,SRR6797128",
                    "coassembly_1"
                ]),
                ""
            ]
        )
        with open(cluster_path) as f:
            self.assertEqual(expected, f.read())

    def test_update_sra_download_real(self):
        output_dir = os.path.join("example", "test_update_sra_download_real")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken update "
            f"--assemble-unmapped "
            f"--forward SRR8334323 SRR8334324 SRR7039260 "
            f"--sra "
            f"--genomes {GENOMES} "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--output {output_dir} "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        sra_1_path = os.path.join(output_dir, "coassemble", "sra", "SRR8334323_1.fastq.gz")
        self.assertTrue(os.path.exists(sra_1_path))
        with gzip.open(sra_1_path) as f:
            file = f.readline().decode()
            self.assertTrue("@SRR8334323.1 HS2:487:H80UEADXX:1:1101:1148:1986/1" in file)
            self.assertTrue("@SRR8334323.2 HS2:487:H80UEADXX:1:1101:1148:1986/2" not in file)

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_2.fastq.gz")))

    @pytest.mark.qsub
    def test_update_aviary_run_real(self):
        output_dir = os.path.join("example", "test_update_aviary_run_real")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken update "
            f"--assemble-unmapped "
            f"--forward SRR8334323 SRR8334324 "
            f"--sra "
            f"--run-aviary "
            f"--aviary-speed fast "
            f"--cores 32 "
            f"--assembly-strategy megahit "
            f"--aviary-gtdbtk-db {GTDBTK_DB} "
            f"--aviary-checkm2-db {CHECKM2_DB} "
            f"--aviary-metabuli-db {METABULI_DB} "
            f"--aviary-extra-binners taxvamb "
            f"--genomes {GENOMES} "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--output {output_dir} "
            f"--snakemake-profile aqua "
            f"--local-cores 12 "
            f"--cluster-retries 1 "
            f"--cluster-submission "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_2.fastq.gz")))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta")))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "checkm_minimal.tsv")))

    @pytest.mark.qsub
    def test_update_aviary_run_real_large_assembly(self):
        output_dir = os.path.join("example", "test_update_aviary_run_real_large_assembly")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken update "
            f"--assemble-unmapped "
            f"--forward SRR5753868 SRR5753874 "
            f"--sra "
            f"--genomes {GENOMES} "
            f"--run-aviary "
            f"--cores 32 "
            f"--aviary-gtdbtk-db {GTDBTK_DB} "
            f"--aviary-checkm2-db {CHECKM2_DB} "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra_large.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--output {output_dir} "
            f"--snakemake-profile aqua "
            f"--local-cores 12 "
            f"--cluster-retries 1 "
            f"--cluster-submission "
        )
        subprocess.run(cmd, shell=True, check=True)


        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR5753868_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR5753868_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR5753874_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR5753874_2.fastq.gz")))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta")))

        config_path = os.path.join(output_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "config.yaml")
        self.assertTrue(os.path.exists(config_path))
        with open(config_path) as f:
            config = YAML().load(f)
        self.assertTrue(config["use_megahit"])

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "checkm_minimal.tsv")))

    def test_update_specific_coassembly_sra(self):
        output_dir = os.path.join("example", "test_update_specific_coassembly_sra")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken update "
            f"--forward SRR8334323 SRR8334324 "
            f"--sra "
            f"--genomes {GENOMES} "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra_two.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--coassemblies coassembly_0 "
            f"--output {output_dir} "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334323_2.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_1.fastq.gz")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334324_2.fastq.gz")))

        self.assertFalse(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334320_1.fastq.gz")))
        self.assertFalse(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334320_2.fastq.gz")))
        self.assertFalse(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334321_1.fastq.gz")))
        self.assertFalse(os.path.exists(os.path.join(output_dir, "coassemble", "sra", "SRR8334321_2.fastq.gz")))

    def test_build_with_downloads(self):
        output_dir = os.path.abspath(os.path.join("example", "test_build_with_downloads"))
        self.setup_output_dir(output_dir)

        with in_tempdir():
            path_to_metapackage = os.path.abspath("metapackage.smpkg")
            path_to_checkm2_db = os.path.abspath("checkm2_db")
            path_to_gtdbtk_db = os.path.abspath("gtdb_release")

            cmd = (
                f"binchicken build "
                f"--singlem-metapackage {path_to_metapackage} "
                # f"--gtdbtk-db {path_to_gtdbtk_db} "
                f"--checkm2-db {path_to_checkm2_db} "
                f"--download-databases "
                f"--cluster-retries 0 "
                f"--output {output_dir} "
            )
            subprocess.run(cmd, shell=True, check=True)

            # Check ENV variables
            cmd = "conda env config vars list"
            output = extern.run(cmd).strip().split("\n")

            self.assertTrue(f"SINGLEM_METAPACKAGE_PATH = {path_to_metapackage}" in output)
            # self.assertTrue(f"GTDBTK_DATA_PATH = {path_to_gtdbtk_db}" in output)
            self.assertTrue(f"CHECKM2DB = {path_to_checkm2_db}" in output)
            self.assertTrue(f"TMPDIR = /tmp" in output)

            # Check databases downloaded
            self.assertTrue(os.path.exists(path_to_metapackage))
            # self.assertTrue(os.path.exists(path_to_gtdbtk_db))
            self.assertTrue(os.path.exists(path_to_checkm2_db))

    def test_single_assembly_provided(self):
        output_dir = os.path.join("example", "test_single_assembly_provided")
        self.setup_output_dir(output_dir)

        cmd = (
            f"binchicken single "
            f"--forward {SAMPLE_READS_FORWARD} "
            f"--reverse {SAMPLE_READS_REVERSE} "
            f"--singlem-metapackage {METAPACKAGE} "
            f"--coassembly-samples sample_1 sample_2 "
            f"--prior-assemblies {PRIOR_ASSEMBLY} "
            f"--output {output_dir} "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        cluster_path = os.path.join(output_dir, "coassemble", "target", "elusive_clusters.tsv")
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

        cmd = (
            f"binchicken single "
            f"--forward {SAMPLE_READS_FORWARD} "
            f"--reverse {SAMPLE_READS_REVERSE} "
            f"--singlem-metapackage {METAPACKAGE} "
            f"--coassembly-samples sample_1 sample_2 "
            f"--prior-assemblies {PRIOR_ASSEMBLY} "
            f"--output {output_dir} "
            f"--run-aviary "
            f"--aviary-gtdbtk-db {GTDBTK_DB} "
            f"--aviary-checkm2-db {CHECKM2_DB} "
        )
        try:
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')
        except subprocess.CalledProcessError as e:
            output = e.stderr.decode('ascii')

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        assembly_1_path = os.path.join(output_dir, "coassemble", "coassemble", "sample_1", "assemble", "assembly", "final_contigs.fasta")
        self.assertTrue(os.path.exists(assembly_1_path))
        with open(assembly_1_path) as f:
            lines = f.readlines()
            self.assertEqual(">NODE_1_length_138944_cov_12.417585\n", lines[0])
            self.assertEqual(9700, len(lines))

        assembly_2_path = os.path.join(output_dir, "coassemble", "coassemble", "sample_2", "assemble", "assembly", "final_contigs.fasta")
        self.assertTrue(os.path.exists(assembly_2_path))
        with open(assembly_2_path) as f:
            lines = f.readlines()
            self.assertEqual(">NODE_1_length_138944_cov_12.417585\n", lines[0])
            self.assertEqual(9693, len(lines))

        self.assertTrue("aviary_assemble" not in output)
        self.assertTrue("prior_assemble" in output)
        self.assertTrue("aviary_recover" in output)
        self.assertTrue("aviary_combine" in output)

        print(output)

    def test_update_assembly_provided(self):
        output_dir = os.path.join("example", "test_update_assembly_provided")
        self.setup_output_dir(output_dir)
        update_dir = os.path.join("example", "test_update_assembly_provided_update")
        self.setup_output_dir(update_dir)

        cmd = (
            f"binchicken coassemble "
            f"--forward {SAMPLE_READS_FORWARD_EMPTY} "
            f"--reverse {SAMPLE_READS_REVERSE_EMPTY} "
            f"--genomes {GENOMES} "
            f"--prodigal-meta "
            f"--singlem-metapackage {METAPACKAGE} "
            f"--output {output_dir} "
        )
        subprocess.run(cmd, shell=True, check=True)

        config_path = os.path.join(output_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        cluster_path = os.path.join(output_dir, "coassemble", "target", "elusive_clusters.tsv")
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
                    "sample_1,sample_2",
                    "2",
                    "3",
                    "8758",
                    "sample_1,sample_2,sample_3",
                    "coassembly_0"
                ]),
                ""
            ]
        )
        with open(cluster_path) as f:
            self.assertEqual(expected, f.read())

        cmd = (
            f"binchicken update "
            f"--coassemble-output {output_dir} "
            f"--coassemblies coassembly_0 "
            f"--prior-assemblies {PRIOR_COASSEMBLY} "
            f"--output {update_dir} "
            f"--run-aviary "
            f"--aviary-gtdbtk-db {GTDBTK_DB} "
            f"--aviary-checkm2-db {CHECKM2_DB} "
        )
        try:
            output_raw = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            output = output_raw.stderr.decode('ascii')
        except subprocess.CalledProcessError as e:
            output = e.stderr.decode('ascii')

        config_path = os.path.join(update_dir, "config.yaml")
        self.assertTrue(os.path.exists(config_path))

        assembly_1_path = os.path.join(update_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta")
        self.assertTrue(os.path.exists(assembly_1_path))
        with open(assembly_1_path) as f:
            lines = f.readlines()
            self.assertEqual(">NODE_1_length_138944_cov_12.417585\n", lines[0])
            self.assertEqual(9700, len(lines))

        self.assertTrue("aviary_assemble" not in output)
        self.assertTrue("prior_assemble" in output)
        self.assertTrue("aviary_recover" in output)
        self.assertTrue("aviary_combine" in output)

        print(output)


if __name__ == '__main__':
    unittest.main()
