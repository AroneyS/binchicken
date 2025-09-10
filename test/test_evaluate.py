#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
from snakemake.io import load_configfile
import polars as pl

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

METAPACKAGE = os.path.join(path_to_data, "singlem_metapackage.smpkg")
MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_UNBINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "unbinned.otu_table.tsv")
MOCK_BINNED = os.path.join(MOCK_COASSEMBLE, "coassemble", "appraise", "binned.otu_table.tsv")
MOCK_TARGETS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "targets.tsv")
MOCK_ELUSIVE_EDGES = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_edges.tsv")
MOCK_ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "coassemble", "target", "elusive_clusters.tsv")
MOCK_SUMMARY = os.path.join(MOCK_COASSEMBLE, "coassemble", "summary.tsv")

MOCK_COASSEMBLIES = ' '.join([os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0")])
MOCK_GENOMES = " ".join([
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_1.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_2.fna"),
    os.path.join(MOCK_COASSEMBLE, "coassemble", "coassemble", "coassembly_0", "recover", "bins", "final_bins", "bin_3.fna"),
])

GENOMES = " ".join([os.path.join(path_to_data, "GB_GCA_013286235.1.fna")])
TWO_GENOMES = " ".join([
    os.path.join(path_to_data, "GB_GCA_013286235.1.fna"),
    os.path.join(path_to_data, "GB_GCA_013286235.2.fna"),
    ])

BIN_INFO_COLUMNS = {
    "Bin Id": str,
    "Marker lineage": str,
    "# genomes": str,
    "# markers": str,
    "# marker sets": str,
    "Completeness (CheckM1)": str,
    "Contamination (CheckM1)": str,
    "Strain heterogeneity": str,
    "Genome size (bp)": str,
    "# ambiguous bases": str,
    "# scaffolds": str,
    "# contigs": str,
    "N50 (scaffolds)": str,
    "N50 (contigs)": str,
    "Mean scaffold length (bp)": str,
    "Mean contig length (bp)": str,
    "Longest scaffold (bp)": str,
    "Longest contig (bp)": str,
    "GC": str,
    "GC std (scaffolds > 1kbp)": str,
    "Coding density": str,
    "Translation table": str,
    "# predicted genes": str,
    "0": str,
    "1": str,
    "2": str,
    "3": str,
    "4": str,
    "5+": str,
    "Completeness (CheckM2)": str,
    "Contamination (CheckM2)": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

ELUSIVE_CLUSTERS_COLUMNS = {
    "samples": str,
    "length": int,
    "total_targets": float,
    "total_size": int,
    "recover_samples": str,
    "coassembly": str,
}

CHECKM2_QUALITY_COLUMNS = {
    "Name": str,
    "Completeness": str,
    "Contamination": str,
    "Completeness_Model_Used": str,
    "Translation_Table_Used": str,
    "Coding_Density": str,
    "Contig_N50": str,
    "Average_Gene_Length": str,
    "Genome_Size": str,
    "GC_Content": str,
    "Total_Coding_Sequences": str,
    "Total_Contigs": str,
    "Max_Contig_Length": str,
    "Additional_Notes": str,
}

def write_string_to_file(string, filename):
    with open(filename, "w") as f:
        f.write("\n".join(string.split(" ")))

class Tests(unittest.TestCase):
    def test_evaluate(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--checkm-version 2 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "statistic",
                        "within",
                        "match",
                        "nonmatch",
                        "total",
                        "match_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "recovery",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_bin_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_unbin_sequences",
                        "recovery",
                        "0",
                        "3",
                        "3",
                        "0.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "novel_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-0_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.18", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_specified_files(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("prodigal_bins" in output)
            self.assertTrue("singlem_pipe_bins" in output)
            self.assertTrue("singlem_summarise_bins" in output)
            self.assertTrue("cluster_original_bins" not in output)
            self.assertTrue("cluster_updated_bins" not in output)
            self.assertTrue("summarise_clusters" not in output)
            self.assertTrue("evaluate" in output)
            self.assertTrue("evaluate_plots" in output)

    def test_evaluate_genome_input(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--new-genomes {MOCK_GENOMES} "
                f"--prodigal-meta "
                f"--coassembly-run coassembly_0 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))
            expected = "\n".join(
                [
                    "\t".join([
                        "coassembly",
                        "statistic",
                        "within",
                        "match",
                        "nonmatch",
                        "total",
                        "match_percent",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "sequences",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "taxonomy",
                        "targets",
                        "1",
                        "1",
                        "2",
                        "50.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "bins",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_bin_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "nontarget_unbin_sequences",
                        "recovery",
                        "0",
                        "3",
                        "3",
                        "0.0",
                    ]),
                    "\t".join([
                        "coassembly_0",
                        "novel_sequences",
                        "recovery",
                        "1",
                        "2",
                        "3",
                        "33.33",
                    ]),
                    ""
                ]
            )
            with open(summary_path) as f:
                self.assertEqual(expected, f.read())

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_1_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_3_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-bin_3_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.14", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_genome_files_input(self):
        with in_tempdir():
            write_string_to_file(MOCK_GENOMES, "genomes")

            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--new-genomes-list genomes "
                f"--prodigal-meta "
                f"--coassembly-run coassembly_0 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("prodigal_bins" in output)
            self.assertTrue("singlem_pipe_bins" in output)
            self.assertTrue("singlem_summarise_bins" in output)
            self.assertTrue("cluster_original_bins" not in output)
            self.assertTrue("cluster_updated_bins" not in output)
            self.assertTrue("summarise_clusters" not in output)
            self.assertTrue("evaluate" in output)
            self.assertTrue("evaluate_plots" in output)

    def test_evaluate_default_config(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
            )
            extern.run(cmd)

            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertEqual(config["max_threads"], 8)
            self.assertEqual(config["checkm_version"], 2)
            self.assertEqual(config["min_completeness"], 70)
            self.assertEqual(config["max_contamination"], 10)

    def test_evaluate_metapackage_env_variable(self):
        with in_tempdir():
            os.environ['SINGLEM_METAPACKAGE_PATH'] = METAPACKAGE
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--output test "
                f"--snakemake-args \"singlem_summarise_bins\" "
            )
            import subprocess
            _ = subprocess.run(cmd, shell=True, check=True, capture_output=True, env=os.environ)

            summarise_path = os.path.join("test", "evaluate", "summarise", "bins_summarised.otu_table.tsv")
            self.assertTrue(os.path.exists(summarise_path))
            with open(summarise_path) as f:
                observed = f.read()

            otu1 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-0_transcripts", "TATCAAGTTCCACAAGAAGTTAGAGGAGAAAGAAGAATCTCGTTAGCTATTAGATGGATT", "1", "1.14", "Root; d__Bacteria"])
            otu2 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TACCAGGTTCCTGTTGAAGTTCGCCCATCACGCCGTTCTGCTTTGGCAATGCGCTGGGTG", "1", "1.14", "Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Polynucleobacter; s__Polynucleobacter_sp018688335"])
            otu3 = "\t".join(["S3.7.ribosomal_protein_S7", "coassembly_0-1_transcripts", "TTCCAGGTGCCTACTGAAGTTCGTCCCGAGCGTAAAATTGCATTGGGTATGAAATGGCTG", "1", "1.18", "Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Mucilaginibacter; s__Mucilaginibacter_sp013286235"])

            self.assertTrue(otu1 in observed)
            self.assertTrue(otu2 in observed)
            self.assertTrue(otu3 in observed)

    def test_evaluate_cluster(self):
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--aviary-outputs {MOCK_COASSEMBLIES} "
                f"--cluster "
                f"--genomes {TWO_GENOMES} "
                f"--checkm-version 2 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            summary_path = os.path.join("test", "evaluate", "evaluate", "summary_stats.tsv")
            self.assertTrue(os.path.exists(summary_path))

            cluster_path = os.path.join("test", "evaluate", "evaluate", "cluster_stats.csv")
            self.assertTrue(os.path.exists(cluster_path))
            expected = "\n".join(
                [
                    ",".join([
                        "original",
                        "2",
                    ]),
                    ",".join([
                        "coassembly_0",
                        "4",
                    ]),
                    ""
                ]
            )
            with open(cluster_path) as f:
                self.assertEqual(expected, f.read())

    def test_evaluate_prefers_recovered_bins(self):
        with in_tempdir():
            # elusive clusters
            elusive_clusters_path = os.path.join("prior_coassemble", "coassemble", "target", "elusive_clusters.tsv")
            os.makedirs(os.path.dirname(elusive_clusters_path), exist_ok=True)
            (
                pl.DataFrame([
                        ["sample_1,sample_2", "2", "2", "2869", "sample_1,sample_2", "coassembly_0"],
                        ["sample_1,sample_3", "2", "2", "2869", "sample_1,sample_3", "coassembly_1"],
                        ],
                    orient="row", schema=ELUSIVE_CLUSTERS_COLUMNS
                    )
                .write_csv(elusive_clusters_path, separator="\t")
            )

            unbinned_path = os.path.join("prior_coassemble", "coassemble", "appraise", "unbinned.otu_table.tsv")
            binned_path = os.path.join("prior_coassemble", "coassemble", "appraise", "binned.otu_table.tsv")
            os.makedirs(os.path.dirname(unbinned_path), exist_ok=True)
            read_size_path = os.path.join("prior_coassemble", "coassemble", "read_size.csv")

            with open(unbinned_path, "w") as f:
                pass
            with open(binned_path, "w") as f:
                pass
            with open(read_size_path, "w") as f:
                pass

            # recovered_bins outputs
            output_bins_path = os.path.join("prior_coassemble", "recovered")
            os.makedirs(os.path.join(output_bins_path, "bins"), exist_ok=True)

            genome1_path = os.path.join(output_bins_path, "bins", "binchicken_coassembly_0.1.fna")
            with open(genome1_path, "w") as f:
                f.write(">contigA\nACGT\n")

            genome2_path = os.path.join(output_bins_path, "bins", "binchicken_coassembly_0.2.fna")
            with open(genome2_path, "w") as f:
                f.write(">contigB\nTGCA\n")

            bin_info = pl.DataFrame([
                    [
                        "binchicken_coassembly_0.1",
                        "p__Bacteroidetes (UID2605)", "350", "314", "208", "86.44", "1.8", "12.5", "1705643", "0", "315", "315", "6824", "6824", "5414", "5414", "25260", "25260", "59.9", "3.41", "90.22", "11",
                        "1669", "37", "269", "8", "0", "0", "0",
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "binchicken_coassembly_0.2",
                        "o__Clostridiales (UID1212)", "172", "257", "149", "44.16", "1.03", "33.33", "980304", "0", "399", "399", "2442", "2442", "2456", "2456", "7735", "7735", "60.2", "2.97", "90.94", "11",
                        "1244", "143", "111", "3", "0", "0", "0",
                        "36.49", "0.24", "Neural Network (Specific Model)", "11", "0.909", "2442", "239.274115755627", "980304", "0.6", "1244", "399", "7735", ""
                    ],
                ], orient="row", schema=BIN_INFO_COLUMNS)
            (
                bin_info
                .write_csv(os.path.join(output_bins_path, "bin_info.tsv"), separator="\t")
            )

            (
                bin_info
                .rename({
                    "Bin Id": "Name",
                    "Completeness (CheckM2)": "Completeness",
                    "Contamination (CheckM2)": "Contamination",
                    })
                .select(CHECKM2_QUALITY_COLUMNS.keys())
                .write_csv(os.path.join(output_bins_path, "quality_report.tsv"), separator="\t")
            )

            # competing coassembly_1 outputs
            source_bins_path = os.path.join("test", "coassemble", "coassemble", "coassembly_1", "recover", "bins")
            os.makedirs(os.path.join(source_bins_path, "final_bins"), exist_ok=True)

            source_genome1_path = os.path.join(source_bins_path, "final_bins", "metabat2_refined_bins.tsv.005.fna")
            with open(source_genome1_path, "w") as f:
                f.write(">contigA\nACGT\n")

            source_genome2_path = os.path.join(source_bins_path, "final_bins", "semibin_refined_bins.tsv.072_sub.fna")
            with open(source_genome2_path, "w") as f:
                f.write(">contigB\nTGCA\n")

            bin_info = pl.DataFrame([
                    [
                        "metabat2_refined_bins.tsv.005",
                        "p__Bacteroidetes (UID2605)", "350", "314", "208", "86.44", "1.8", "12.5", "1705643", "0", "315", "315", "6824", "6824", "5414", "5414", "25260", "25260", "59.9", "3.41", "90.22", "11",
                        "1669", "37", "269", "8", "0", "0", "0",
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "semibin_refined_bins.tsv.072_sub",
                        "o__Clostridiales (UID1212)", "172", "257", "149", "44.16", "1.03", "33.33", "980304", "0", "399", "399", "2442", "2442", "2456", "2456", "7735", "7735", "60.2", "2.97", "90.94", "11",
                        "1244", "143", "111", "3", "0", "0", "0",
                        "36.49", "0.24", "Neural Network (Specific Model)", "11", "0.909", "2442", "239.274115755627", "980304", "0.6", "1244", "399", "7735", ""
                    ],
                ], orient="row", schema=BIN_INFO_COLUMNS)
            (
                bin_info
                .write_csv(os.path.join(source_bins_path, "bin_info.tsv"), separator="\t")
            )

            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output prior_coassemble "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            # Assert: pipeline stages that consume recovered bins are present
            self.assertTrue("singlem_summarise_bins" in output)
            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertIn("recovered_bins", config)
            self.assertEqual(sorted(config["recovered_bins"].keys()), ["binchicken_coassembly_0.1","binchicken_coassembly_0.2"])


if __name__ == '__main__':
    unittest.main()
