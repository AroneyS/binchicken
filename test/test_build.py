#!/usr/bin/env python3

import unittest
import os
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
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])

class Tests(unittest.TestCase):
    def test_build(self):
        with in_tempdir():
            #path_to_conda = os.path.abspath(".conda")
            cmd = (
                f"binchicken build "
                f"--conda-prefix {path_to_conda} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--gtdbtk-db gtdb_release "
                f"--checkm2-db CheckM2_database "
                f"--set-tmp-dir /tmp "
                f"--skip-aviary-envs "
            )
            extern.run(cmd)

            # Check ENV variables
            cmd = "conda env config vars list"
            output = extern.run(cmd).strip().split("\n")

            self.assertTrue(f"SNAKEMAKE_CONDA_PREFIX = {path_to_conda}" in output)
            self.assertTrue(f"CONDA_ENV_PATH = {path_to_conda}" in output)
            self.assertTrue(f"SINGLEM_METAPACKAGE_PATH = {METAPACKAGE}" in output)
            self.assertTrue(f"GTDBTK_DATA_PATH = gtdb_release" in output)
            self.assertTrue(f"CHECKM2DB = CheckM2_database" in output)
            self.assertTrue(f"TMPDIR = /tmp" in output)

            # Dryrun coassemble
            cmd = (
                f"binchicken coassemble "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
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
                f"--singlem-metapackage {METAPACKAGE} "
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
                f"--aviary-checkm2-db CheckM2_database "
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
                f"--singlem-metapackage {METAPACKAGE} "
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
