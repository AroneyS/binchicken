#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import pandas as pd
import snakemake as smk
from ruamel.yaml import YAML

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')
path_to_ibis = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "ibis")
path_to_smk = os.path.join(path_to_ibis, "workflow", "evaluate.smk")
path_to_config = os.path.join(path_to_ibis, "config", "template_evaluate.yaml")

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
TARGETS = os.path.join(MOCK_COASSEMBLE, "target", "targets.tsv")
BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
ELUSIVE_EDGES = os.path.join(MOCK_COASSEMBLE, "target", "elusive_edges.tsv")
ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters.tsv")

MATCH_COLUMNS = ["coassembly", "gene", "sequence", "genome", "target", "found_in", "taxonomy"]
SUM_STATS_COLUMNS = ["coassembly", "statistic", "within", "match", "nonmatch", "total", "match_percent"]
COASSEMBLE_SUM_COLUMNS = ["coassembly", "samples", "length", "total_weight", "total_targets", "total_size", "unmapped_size"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def run_snakemake(self, matched_hits, novel_hits, coassemble_summary, recovered_bins, cluster_summary=None, cluster_ani=0.86):
        # copy inputs into temp folders
        os.makedirs(os.path.join(os.getcwd(), "evaluate", "evaluate"))
        matched_hits.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "matched_hits.tsv"), sep="\t")
        novel_hits.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "novel_hits.tsv"), sep="\t")
        if cluster_summary:
            cluster_summary.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "cluster_stats.tsv"), sep="\t")

        coassemble_summary_path = os.path.join(os.getcwd(), "coassemble_summary.tsv")
        coassemble_summary.to_csv(coassemble_summary_path, sep="\t")

        # Load config
        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False
        with open(path_to_config) as f:
            test_config = yaml.load(f)
        test_config["test"] = True
        test_config["recovered_bins"] = recovered_bins
        test_config["singlem_metapackage"] = METAPACKAGE
        test_config["targets"] = TARGETS
        test_config["binned"] = BINNED
        test_config["elusive_edges"] = ELUSIVE_EDGES
        test_config["elusive_clusters"] = ELUSIVE_CLUSTERS
        test_config["coassemble_summary"] = coassemble_summary_path
        test_config["cluster"] = cluster_ani if cluster_summary else False

        # Run snakemake
        updated_files = []
        smk.snakemake(
            path_to_smk, config=test_config,
            targets=["evaluate_plots"],
            updated_files=updated_files,
            use_conda=True, conda_prefix=path_to_conda
            )

        # load outputs from temp folders and return
        outputs = {os.path.basename(f): pd.read_csv(f, sep="\t") for f in updated_files if f.endswith(".tsv")}
        return outputs

    def test_evaluate_Rscript(self):
        with in_tempdir():
            matched_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
                ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
                ["coassembly_0", "S3.1", "ZZZ", None, "99", None, "Root"],
                ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
                ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
                ["coassembly_1", "S3.1", "YYY", None, "98", None, "Root"],
            ], columns=MATCH_COLUMNS)
            novel_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
                ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
            ], columns=MATCH_COLUMNS)
            coassemble_summary = pd.DataFrame([
                ["coassembly_0", "sample_1,sample_2,sample_3", 3, 2, 2, 3000, 300],
                ["coassembly_1", "sample_1,sample_2", 2, 2, 2, 2000, 200],
            ], columns=COASSEMBLE_SUM_COLUMNS)
            recovered_bins = {
                "coassembly_0-genome_0": "genome_0.fna",
                "coassembly_0-genome_1": "genome_1.fna",
                "coassembly_1-genome_0": "genome_0.fna",
                "coassembly_1-genome_1": "genome_1.fna",
                "coassembly_1-genome_2": "genome_2.fna",
                }

            expected_sum_stats = pd.DataFrame([
                ["coassembly_0", "sequences", "targets", 1, 1, 2, 50],
                ["coassembly_0", "taxonomy", "targets", 1, 1, 2, 50],
                ["coassembly_0", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_0", "novel_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_0", "bins", "recovery", 1, 1, 2, 50],
                ["coassembly_1", "sequences", "targets", 1, 1, 2, 50],
                ["coassembly_1", "taxonomy", "targets", 1, 1, 2, 50],
                ["coassembly_1", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_1", "bins", "recovery", 1, 2, 3, 33.33],
            ], columns=SUM_STATS_COLUMNS)
            observed = self.run_snakemake(matched_hits, novel_hits, coassemble_summary, recovered_bins)
            observed_sum_stats = observed["summary_stats.tsv"]
            self.assertDataFrameEqual(expected_sum_stats, observed_sum_stats)

    def test_evaluate_Rscript_duplicates(self):
        with in_tempdir():
            matched_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAA", "genome_0_transcripts", "10", None, "Root"],
                ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
                ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
                ["coassembly_0", "S3.1", "ZZZ", None, "99", None, "Root"],
                ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
                ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
                ["coassembly_1", "S3.1", "YYY", None, "98", None, "Root"],
            ], columns=MATCH_COLUMNS)
            novel_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
                ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
            ], columns=MATCH_COLUMNS)
            coassemble_summary = pd.DataFrame([
                ["coassembly_0", "sample_1,sample_2,sample_3", 3, 2, 2, 3000, 300],
                ["coassembly_1", "sample_1,sample_2", 2, 2, 2, 2000, 200],
            ], columns=COASSEMBLE_SUM_COLUMNS)
            recovered_bins = {
                "coassembly_0-genome_0": "genome_0.fna",
                "coassembly_0-genome_1": "genome_1.fna",
                "coassembly_1-genome_0": "genome_0.fna",
                "coassembly_1-genome_1": "genome_1.fna",
                "coassembly_1-genome_2": "genome_2.fna",
                }

            expected_sum_stats = pd.DataFrame([
                ["coassembly_0", "sequences", "targets", 1, 1, 2, 50],
                ["coassembly_0", "taxonomy", "targets", 1, 1, 2, 50],
                ["coassembly_0", "nontarget_sequences", "recovery", 1, 3, 4, 25],
                ["coassembly_0", "novel_sequences", "recovery", 1, 3, 4, 25],
                ["coassembly_0", "bins", "recovery", 2, 0, 2, 100],
                ["coassembly_1", "sequences", "targets", 1, 1, 2, 50],
                ["coassembly_1", "taxonomy", "targets", 1, 1, 2, 50],
                ["coassembly_1", "nontarget_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_1", "novel_sequences", "recovery", 1, 2, 3, 33.33],
                ["coassembly_1", "bins", "recovery", 1, 2, 3, 33.33],
            ], columns=SUM_STATS_COLUMNS)
            observed = self.run_snakemake(matched_hits, novel_hits, coassemble_summary, recovered_bins)
            observed_sum_stats = observed["summary_stats.tsv"]
            self.assertDataFrameEqual(expected_sum_stats, observed_sum_stats)

if __name__ == '__main__':
    unittest.main()
