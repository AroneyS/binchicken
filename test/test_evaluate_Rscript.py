#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import pandas as pd
import snakemake as smk
from ruamel.yaml import YAML
import extern

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

MATCH_COLUMNS = ["coassembly", "gene", "sequence", "genome", "target", "found_in", "taxonomy"]
SUM_STATS_COLUMNS = ["coassembly", "statistic", "within", "match", "nonmatch", "total", "match_percent"]

class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        pd.testing.assert_frame_equal(a, b)

    def run_snakemake(self, matched_hits, novel_hits, coassemble_summary=None, cluster_summary=None, cluster_ani=0.86):
        # copy inputs into temp folders
        os.makedirs(os.path.join(os.getcwd(), "evaluate", "evaluate"))
        matched_hits.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "matched_hits.tsv"), sep="\t")
        novel_hits.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "novel_hits.tsv"), sep="\t")
        if cluster_summary:
            cluster_summary.to_csv(os.path.join(os.getcwd(), "evaluate", "evaluate", "cluster_stats.tsv"), sep="\t")

        # Load config
        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False
        with open(path_to_config) as f:
            test_config = yaml.load(f)
        test_config["test"] = True
        test_config["recovered_bins"] = {
                "coassembly_0-0": '/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_1.fna',
                "coassembly_0-1": "/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_3.fna",
                }
        test_config["singlem_metapackage"] = METAPACKAGE
        test_config["targets"] = TARGETS
        test_config["binned"] = BINNED
        test_config["elusive_edges"] = ELUSIVE_EDGES
        test_config["elusive_clusters"] = ELUSIVE_CLUSTERS
        test_config["coassemble_summary"] = coassemble_summary if coassemble_summary else COASSEMBLE_SUMMARY
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

    def run_rscript(self, matched_hits, novel_hits, coassemble_summary=None, cluster_summary=None):
        matched_hits_path = os.path.join(os.getcwd(), "matched_hits.tsv")
        matched_hits.to_csv(matched_hits_path, sep="\t")

        novel_hits_path = os.path.join(os.getcwd(), "novel_hits.tsv")
        novel_hits.to_csv(novel_hits_path, sep="\t")

        if cluster_summary:
            cluster_summary_path = os.path.join(os.getcwd(), "cluster_stats.tsv")
            cluster_summary.to_csv(cluster_summary_path, sep="\t")
        else:
            cluster_summary_path = "NA"

        if coassemble_summary:
            coassemble_summary_path = os.path.join(os.getcwd(), "coassemble_summary.tsv")
            coassemble_summary.to_csv(coassemble_summary_path, sep="\t")
        else:
            coassemble_summary_path = COASSEMBLE_SUMMARY

        recovered_bins = {
                "coassembly_0-0": '/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_1.fna',
                "coassembly_0-1": "/mnt/hpccs01/scratch/microbiome/aroneys/src/cockatoo_dev/test/data/mock_coassemble/coassemble/coassembly_0/recover/bins/final_bins/bin_3.fna",
                }

        main_dir = os.path.join(os.getcwd(), "plots")
        summary_stats_path = os.path.join(os.getcwd(), "summary_stats.tsv")
        summary_table_path = os.path.join(os.getcwd(), "summary_table.png")

        cmd = (
            f"conda run -p {path_to_conda}/a29444b44cc7faccbc9bb084f77b29b0_ "
            f"Rscript "
            f"{path_to_cockatoo}/workflow/scripts/evaluate.R "
            f"matched_hits_path={matched_hits_path} "
            f"novel_hits_path={novel_hits_path} "
            f"cluster_summary_path={cluster_summary_path} "
            f"coassemble_summary_path={coassemble_summary_path} "
            f"main_dir={main_dir} "
            f"summary_stats_path={summary_stats_path} "
            f"summary_table_path={summary_table_path} "
        )
            #f"recovered_bins={recovered_bins} "
        extern.run(cmd)


    def test_snakemake(self):
        with in_tempdir():
            matched_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAA", "genome_1_transcripts", "10", None, "Root"],
                ["coassembly_0", "S3.1", "BBB", "genome_1_transcripts", None, "oldgenome_1", "Root"],
                ["coassembly_1", "S3.1", "CCC", "genome_1_transcripts", "11", None, "Root"],
                ["coassembly_1", "S3.1", "DDD", "genome_1_transcripts", None, "oldgenome_2", "Root"],
            ], columns=MATCH_COLUMNS)
            novel_hits = pd.DataFrame([
                ["coassembly_0", "S3.1", "AAB", "genome_1_transcripts", None, None, "Root"],
                ["coassembly_1", "S3.1", "CCD", "genome_1_transcripts", None, None, "Root"],
            ], columns=MATCH_COLUMNS)

            expected_sum_stats = pd.DataFrame([
                ["coassembly_0", "sequences", "targets", "1", "1", "2", "50"],
                ["coassembly_0", "taxonomy", "targets", "1", "1", "2", "50"],
                ["coassembly_0", "nontarget_sequences", "recovery", "1", "2", "3", "33.33"],
                ["coassembly_0", "novel_sequences", "recovery", "1", "2", "3", "33.33"],
                ["coassembly_0", "bins", "recovery", "1", "1", "2", "50"],
                ["coassembly_1", "sequences", "targets", "1", "1", "2", "50"],
                ["coassembly_1", "taxonomy", "targets", "1", "1", "2", "50"],
                ["coassembly_1", "nontarget_sequences", "recovery", "1", "2", "3", "33.33"],
                ["coassembly_1", "novel_sequences", "recovery", "1", "2", "3", "33.33"],
                ["coassembly_1", "bins", "recovery", "1", "1", "2", "50"],
            ], columns=SUM_STATS_COLUMNS)
            observed = self.run_snakemake(matched_hits, novel_hits)
            print(observed["summary_stats.tsv"])
            observed_sum_stats = observed["summary_stats.tsv"]
            #self.assertDataFrameEqual(expected_sum_stats, observed_sum_stats)


if __name__ == '__main__':
    unittest.main()
