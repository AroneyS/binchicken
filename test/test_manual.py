#!/usr/bin/env python3

import unittest
import os
import shutil
import gzip
import subprocess
from ruamel.yaml import YAML
import extern
from bird_tool_utils import in_tempdir

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

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble", "coassemble")
APPRAISE_BINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "binned.otu_table.tsv")
APPRAISE_UNBINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned.otu_table.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_clusters.tsv")
ELUSIVE_CLUSTERS_TWO = os.path.join(MOCK_COASSEMBLE, "coassemble", 'target', 'elusive_clusters_two.tsv')

MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])

class Tests(unittest.TestCase):
    def setup_output_dir(self, output_dir):
        try:
            shutil.rmtree(output_dir)
        except FileNotFoundError:
            pass
        os.makedirs(output_dir)

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
            f"--conda-prefix {path_to_conda} "
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
            f"--aviary-gtdbtk-db /work/microbiome/db/gtdb/gtdb_release207_v2 "
            f"--aviary-checkm2-db /work/microbiome/db/CheckM2_database "
            f"--genomes {GENOMES} "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--output {output_dir} "
            f"--conda-prefix {path_to_conda} "
            f"--snakemake-profile mqsub "
            f"--local-cores 5 "
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
            f"--aviary-gtdbtk-db /work/microbiome/db/gtdb/gtdb_release207_v2 "
            f"--aviary-checkm2-db /work/microbiome/db/CheckM2_database "
            f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'unbinned_sra.otu_table.tsv')} "
            f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'appraise', 'binned_sra.otu_table.tsv')} "
            f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'target', 'targets.tsv')} "
            f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_edges.tsv')} "
            f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_sra_large.tsv')} "
            f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'summary.tsv')} "
            f"--output {output_dir} "
            f"--conda-prefix {path_to_conda} "
            f"--snakemake-profile mqsub "
            f"--local-cores 5 "
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
            f"--conda-prefix {path_to_conda} "
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
        with in_tempdir():
            # path_to_conda = os.path.abspath(".conda")
            path_to_metapackage = "metapackage.smpkg"
            path_to_checkm2_db = "checkm2_db"
            path_to_gtdbtk_db = "gtdb_release"

            cmd = (
                f"binchicken build "
                f"--conda-prefix {path_to_conda} "
                f"--singlem-metapackage {path_to_metapackage} "
                # f"--gtdbtk-db {path_to_gtdbtk_db} "
                f"--checkm2-db {path_to_checkm2_db} "
                f"--download-databases "
            )
            extern.run(cmd)

            # Check ENV variables
            cmd = "conda env config vars list"
            output = extern.run(cmd).strip().split("\n")

            self.assertTrue(f"SNAKEMAKE_CONDA_PREFIX = {path_to_conda}" in output)
            self.assertTrue(f"CONDA_ENV_PATH = {path_to_conda}" in output)
            self.assertTrue(f"SINGLEM_METAPACKAGE_PATH = {path_to_metapackage}" in output)
            # self.assertTrue(f"GTDBTK_DATA_PATH = {path_to_gtdbtk_db}" in output)
            self.assertTrue(f"CHECKM2DB = {path_to_checkm2_db}" in output)
            self.assertTrue(f"TMPDIR = /tmp" in output)

            # Check databases downloaded
            import pdb; pdb.set_trace()

            # Dryrun coassemble
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {path_to_metapackage} "
                f"--output test_coassemble "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            output = extern.run(cmd)
            self.assertFalse("binchicken/workflow/env/singlem.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/coverm.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/r.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/kingfisher.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/prodigal.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/aviary.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/fastp.yml will be created." in output)

            # Dryrun evaluate
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {path_to_metapackage} "
                f"--output test_evaluate "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            output = extern.run(cmd)
            self.assertFalse("binchicken/workflow/env/singlem.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/coverm.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/r.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/kingfisher.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/prodigal.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/aviary.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/fastp.yml will be created." in output)

            # Dryrun update
            cmd = (
                f"binchicken update "
                f"--assemble-unmapped "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--run-aviary "
                f"--aviary-gtdbtk-db gtdb_release "
                f"--aviary-checkm2-db {path_to_checkm2_db} "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'appraise', 'unbinned_sra.otu_table.tsv')} "
                f"--coassemble-binned {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'appraise', 'binned_sra.otu_table.tsv')} "
                f"--coassemble-targets {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'target', 'targets.tsv')} "
                f"--coassemble-elusive-edges {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'target', 'elusive_edges.tsv')} "
                f"--coassemble-elusive-clusters {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'target', 'elusive_clusters_sra.tsv')} "
                f"--coassemble-summary {os.path.join(MOCK_COASSEMBLE, 'coassemble', 'summary.tsv')} "
                f"--output test_update "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            output = extern.run(cmd)
            self.assertFalse("binchicken/workflow/env/singlem.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/coverm.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/r.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/kingfisher.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/prodigal.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/aviary.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/fastp.yml will be created." in output)

            # Dryrun iterate
            cmd = (
                f"binchicken iterate "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {path_to_metapackage} "
                f"--output test_iterate "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
            )
            output = extern.run(cmd)
            self.assertFalse("binchicken/workflow/env/singlem.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/coverm.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/r.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/kingfisher.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/prodigal.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/aviary.yml will be created." in output)
            self.assertFalse("binchicken/workflow/env/fastp.yml will be created." in output)


if __name__ == '__main__':
    unittest.main()
