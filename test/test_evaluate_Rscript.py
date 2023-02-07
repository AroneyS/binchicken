#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import pandas as pd
import snakemake as smk
from ruamel.yaml import YAML

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')
path_to_cockatoo = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "cockatoo")
path_to_smk = os.path.join(path_to_cockatoo, "workflow", "evaluate.smk")
path_to_config = os.path.join(path_to_cockatoo, "config", "template_evaluate.yaml")

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])
METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
TARGETS = os.path.join(MOCK_COASSEMBLE, "target", "targets.tsv")
BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
ELUSIVE_EDGES = os.path.join(MOCK_COASSEMBLE, "target", "elusive_edges.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters.tsv")
COASSEMBLE_SUMMARY = os.path.join(MOCK_COASSEMBLE, "summary.tsv")

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def run_snakemake(self):
        # copy inputs into temp folders

        # Load config
        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False
        with open(path_to_config) as f:
            test_config = yaml.load(f)
        test_config["recovered_bins"] = {
                "coassembly_0-0": '/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_1.fna',
                "coassembly_0-1": "/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_3.fna",
                }
        test_config["singlem_metapackage"] = METAPACKAGE
        test_config["targets"] = TARGETS
        test_config["binned"] = BINNED
        test_config["elusive_edges"] = ELUSIVE_EDGES
        test_config["elusive_clusters"] = ELUSIVE_CLUSTERS
        test_config["coassemble_summary"] = COASSEMBLE_SUMMARY

        # Run snakemake
        updated_files = []
        smk.snakemake(path_to_smk, config=test_config, targets=["evaluate_plots"], updated_files=updated_files, use_conda=True, conda_prefix=path_to_conda)

        # load outputs from temp folders and return

    def test_snakemake(self):
        with in_tempdir():
            self.run_snakemake()


if __name__ == '__main__':
    unittest.main()
