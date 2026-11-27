#!/usr/bin/env python3

import unittest
import os
from bird_tool_utils import in_tempdir
import extern
from snakemake.common.configfile import load_configfile
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
                f"--snakemake-args \" --quiet rules\" "
            )
            output = extern.run(cmd)

            self.assertTrue("\nprodigal_bins" in output)
            self.assertTrue("\nsinglem_pipe_bins" in output)
            self.assertTrue("\nsinglem_summarise_bins" in output)
            self.assertTrue("\ncluster_original_bins" not in output)
            self.assertTrue("\ncluster_updated_bins" not in output)
            self.assertTrue("\nsummarise_clusters" not in output)
            self.assertTrue("\nevaluate" in output)
            self.assertTrue("\nevaluate_plots" in output)

    def test_evaluate_genome_input(self):
        # --new-genomes bins are named <coassembly-run>-<genome> and all evaluated (no quality filtering).
        # The evaluation itself is covered by test_evaluate, so this checks the bins passed to it
        with in_tempdir():
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--new-genomes {MOCK_GENOMES} "
                f"--prodigal-meta "
                f"--coassembly-run coassembly_0 "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            expected = {
                f"coassembly_0-{os.path.splitext(os.path.basename(g))[0]}": os.path.abspath(g)
                for g in MOCK_GENOMES.split(" ")
                }
            self.assertEqual(expected, config["recovered_bins"])

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
                f"--snakemake-args \" --quiet rules\" "
            )
            output = extern.run(cmd)

            self.assertTrue("\nprodigal_bins" in output)
            self.assertTrue("\nsinglem_pipe_bins" in output)
            self.assertTrue("\nsinglem_summarise_bins" in output)
            self.assertTrue("\ncluster_original_bins" not in output)
            self.assertTrue("\ncluster_updated_bins" not in output)
            self.assertTrue("\nsummarise_clusters" not in output)
            self.assertTrue("\nevaluate" in output)
            self.assertTrue("\nevaluate_plots" in output)

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
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "binchicken_coassembly_0.2",
                        "36.49", "0.24", "Neural Network (Specific Model)", "11", "0.909", "2442", "239.274115755627", "980304", "0.6", "1244", "399", "7735", ""
                    ],
                ], orient="row", schema=BIN_INFO_COLUMNS)
            (
                bin_info
                .write_csv(os.path.join(output_bins_path, "bin_info.tsv"), separator="\t")
            )

            (
                bin_info
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
                        "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7531455961654", "1705643", "0.6", "1669", "315", "25260", ""
                    ],
                    [
                        "semibin_refined_bins.tsv.072_sub",
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
                f"--snakemake-args \" --quiet rules\" "
            )
            output = extern.run(cmd)

            # Assert: pipeline stages that consume recovered bins are present
            self.assertTrue("\nsinglem_summarise_bins" in output)
            config = load_configfile(os.path.join("test", "config.yaml"))
            self.assertIn("recovered_bins", config)
            # binchicken_coassembly_0.2 (36.49% complete) fails the default --min-completeness 70
            self.assertEqual(sorted(config["recovered_bins"].keys()), ["binchicken_coassembly_0.1"])


    def _make_recovered_coassemble(self, coassemblies):
        """Mock a prior coassemble run with one passing recovered bin per coassembly."""
        target_dir = os.path.join("prior_coassemble", "coassemble", "target")
        appraise_dir = os.path.join("prior_coassemble", "coassemble", "appraise")
        os.makedirs(target_dir, exist_ok=True)
        os.makedirs(appraise_dir, exist_ok=True)
        (
            pl.DataFrame([
                    ["sample_1,sample_2", "2", "2", "2869", "sample_1,sample_2", c] for c in coassemblies
                    ],
                orient="row", schema=ELUSIVE_CLUSTERS_COLUMNS
                )
            .write_csv(os.path.join(target_dir, "elusive_clusters.tsv"), separator="\t")
        )
        for path in [os.path.join(appraise_dir, "unbinned.otu_table.tsv"), os.path.join(appraise_dir, "binned.otu_table.tsv")]:
            with open(path, "w") as f:
                pass

        recovered_dir = os.path.join("prior_coassemble", "recovered")
        os.makedirs(os.path.join(recovered_dir, "bins"), exist_ok=True)
        rows = []
        for c in coassemblies:
            name = f"binchicken_{c}.1"
            with open(os.path.join(recovered_dir, "bins", name + ".fna"), "w") as f:
                f.write(">contigA\nACGT\n")
            rows.append([name, "83.86", "1.14", "Neural Network (Specific Model)", "11", "0.902", "6824", "307.7", "1705643", "0.6", "1669", "315", "25260", ""])
        (
            pl.DataFrame(rows, orient="row", schema=CHECKM2_QUALITY_COLUMNS)
            .write_csv(os.path.join(recovered_dir, "quality_report.tsv"), separator="\t")
        )
        return os.path.abspath(os.path.join(recovered_dir, "bins"))

    def _cluster_commands(self, output):
        """Map cluster definition name -> planned coverm cluster command line."""
        commands = {}
        for line in output.splitlines():
            if "coverm cluster" not in line:
                continue
            name = line.split("--output-cluster-definition ")[1].split(" ")[0]
            commands[os.path.basename(name)] = line
        return commands

    def test_evaluate_cluster_coassembly_name_prefix(self):
        # coassembly_1 must not pick up bins from coassembly_10 (substring match)
        with in_tempdir():
            bins_dir = self._make_recovered_coassemble(["coassembly_1", "coassembly_10"])
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output prior_coassemble "
                f"--cluster "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --printshellcmds\" "
            )
            commands = self._cluster_commands(extern.run(cmd))

            self.assertEqual(["coassembly_1.txt", "coassembly_10.txt", "original.txt"], sorted(commands))
            self.assertIn(os.path.join(bins_dir, "binchicken_coassembly_1.1.fna"), commands["coassembly_1.txt"])
            self.assertNotIn(os.path.join(bins_dir, "binchicken_coassembly_10.1.fna"), commands["coassembly_1.txt"])
            self.assertIn(os.path.join(bins_dir, "binchicken_coassembly_10.1.fna"), commands["coassembly_10.txt"])
            self.assertNotIn(os.path.join(bins_dir, "binchicken_coassembly_1.1.fna"), commands["coassembly_10.txt"])

    def test_evaluate_cluster_single_assembly_names(self):
        # Single-assembly runs name "coassemblies" after samples, with no coassembly_N in bin names
        with in_tempdir():
            bins_dir = self._make_recovered_coassemble(["sample_1", "sample_2"])
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output prior_coassemble "
                f"--cluster "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --printshellcmds\" "
            )
            commands = self._cluster_commands(extern.run(cmd))

            self.assertEqual(["original.txt", "sample_1.txt", "sample_2.txt"], sorted(commands))
            self.assertIn(os.path.join(bins_dir, "binchicken_sample_1.1.fna"), commands["sample_1.txt"])
            self.assertNotIn(os.path.join(bins_dir, "binchicken_sample_2.1.fna"), commands["sample_1.txt"])


    def test_evaluate_cluster_iteration_like_sample_name(self):
        # A sample named it2_sample must not be read as iteration 2 of "sample"
        with in_tempdir():
            bins_dir = self._make_recovered_coassemble(["it2_sample", "sample"])
            cmd = (
                f"binchicken evaluate "
                f"--coassemble-output prior_coassemble "
                f"--cluster "
                f"--genomes {GENOMES} "
                f"--singlem-metapackage {METAPACKAGE} "
                f"--output test "
                f"--dryrun "
                f"--snakemake-args \" --printshellcmds\" "
            )
            commands = self._cluster_commands(extern.run(cmd))

            self.assertEqual(["it2_sample.txt", "original.txt", "sample.txt"], sorted(commands))
            self.assertIn(os.path.join(bins_dir, "binchicken_it2_sample.1.fna"), commands["it2_sample.txt"])
            self.assertNotIn(os.path.join(bins_dir, "binchicken_it2_sample.1.fna"), commands["sample.txt"])

if __name__ == '__main__':
    unittest.main()
