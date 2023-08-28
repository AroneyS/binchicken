#!/usr/bin/env python3

import unittest
import os
import gzip
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

MOCK_COASSEMBLE = os.path.join(path_to_data, "mock_coassemble")
MOCK_UNBINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "unbinned.otu_table.tsv")
MOCK_BINNED = os.path.join(MOCK_COASSEMBLE, "appraise", "binned.otu_table.tsv")
MOCK_TARGETS = os.path.join(MOCK_COASSEMBLE, "target", "targets.tsv")
MOCK_ELUSIVE_EDGES = os.path.join(MOCK_COASSEMBLE, "target", "elusive_edges.tsv")
MOCK_ELUSIVE_CLUSTERS = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters.tsv")
MOCK_ELUSIVE_CLUSTERS_TWO = os.path.join(MOCK_COASSEMBLE, 'target', 'elusive_clusters_two.tsv')
MOCK_SUMMARY = os.path.join(MOCK_COASSEMBLE, "summary.tsv")

MOCK_UNBINNED_SRA = os.path.join(MOCK_COASSEMBLE, "appraise", "unbinned_sra.otu_table.tsv")
MOCK_BINNED_SRA = os.path.join(MOCK_COASSEMBLE, "appraise", "binned_sra.otu_table.tsv")
MOCK_ELUSIVE_CLUSTERS_SRA = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters_sra.tsv")
MOCK_ELUSIVE_CLUSTERS_SRA_MOCK = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters_sra_mock.tsv")
MOCK_ELUSIVE_CLUSTERS_SRA_MOCK2 = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters_sra_mock2.tsv")
MOCK_ELUSIVE_CLUSTERS_SRA_MOCK3 = os.path.join(MOCK_COASSEMBLE, "target", "elusive_clusters_sra_mock3.tsv")

class Tests(unittest.TestCase):
    def test_update(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            bins_reference_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertFalse(os.path.exists(bins_reference_path))

            output_bam_files = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.bam")
            self.assertFalse(os.path.exists(output_bam_files))

            coverm_working_dir = os.path.join("test", "coassemble", "mapping", "sample_1_coverm")
            self.assertFalse(os.path.exists(coverm_working_dir))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:112:HMNM5DSXX:4:9999:19126:17300" not in file)

            coassemble_path = os.path.join("test", "coassemble", "commands", "coassemble_commands.sh")
            self.assertTrue(os.path.exists(coassemble_path))
            test_dir = os.path.abspath("test")
            expected = "\n".join(
                [
                    " ".join([
                        "aviary assemble -1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble"),
                        "-n 64 -t 64 -m 500 &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_assemble.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(coassemble_path) as f:
                self.assertEqual(expected, f.read())

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            expected = "\n".join(
                [
                    " ".join([
                        "aviary recover --assembly", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "assemble", "assembly", "final_contigs.fasta"),
                        "-1",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.1.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.1.fq.gz"),
                        "-2",
                        os.path.join(test_dir, "coassemble", "mapping", "sample_1_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_2_unmapped.2.fq.gz"),
                        os.path.join(test_dir, "coassemble", "mapping", "sample_3_unmapped.2.fq.gz"),
                        "--output", os.path.join(test_dir, "coassemble", "coassemble", "coassembly_0", "recover"),
                        "--workflow recover_mags_no_singlem --skip-binners maxbin concoct rosella --skip-abundances --refinery-max-iterations 0 "
                        "-n 32 -t 32 -m 250 &>",
                        os.path.join(test_dir, "coassemble", "coassemble", "logs", "coassembly_0_recover.log"),
                        ""
                    ]),
                    ""
                ]
            )
            with open(recover_path) as f:
                self.assertEqual(expected, f.read())

    def test_update_specified_files(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("target_elusive" not in output)
            self.assertTrue("cluster_graph" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" in output)

    def test_update_read_identity(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--unmapping-max-identity 99 "
                f"--unmapping-max-alignment 90 "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            bins_reference_path = os.path.join("test", "coassemble", "mapping", "sample_1_reference.fna")
            self.assertFalse(os.path.exists(bins_reference_path))

            output_bam_files = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.bam")
            self.assertFalse(os.path.exists(output_bam_files))

            coverm_working_dir = os.path.join("test", "coassemble", "mapping", "sample_1_coverm")
            self.assertFalse(os.path.exists(coverm_working_dir))

            unmapped_sample_1_path = os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")
            self.assertTrue(os.path.exists(unmapped_sample_1_path))
            with gzip.open(unmapped_sample_1_path) as f:
                file = f.read().decode()
                self.assertTrue("@A00178:112:HMNM5DSXX:4:1622:16405:19194" in file)
                self.assertTrue("@A00178:118:HTHTVDSXX:1:1249:16740:14105" in file)

    def test_update_sra_download(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output_comb = extern.run(cmd)

            output_sra = output_comb.split("Building DAG of jobs...")[1]
            self.assertTrue("download_sra" in output_sra)
            self.assertTrue("sra_qc" not in output_sra)
            self.assertTrue("compile_sra_qc" not in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output_sra = output_comb.split("Building DAG of jobs...")[2]
            self.assertTrue("download_sra" not in output_sra)
            self.assertTrue("sra_qc" in output_sra)
            self.assertTrue("compile_sra_qc" in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output = output_comb.split("Building DAG of jobs...")[3]
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("target_elusive" not in output)
            self.assertTrue("cluster_graph" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" in output)

    def test_update_sra_download_mock(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--forward SRR3309137 SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED_SRA} "
                f"--coassemble-binned {MOCK_BINNED_SRA} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS_SRA_MOCK} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" --config mock_sra=True\" "
            )
            extern.run(cmd)

            sra_f0_path = os.path.join("test", "coassemble", "sra_qc", "SRR3309137_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f0_path))
            with gzip.open(sra_f0_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" in file)
                self.assertTrue("@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" not in file)

            sra_r0_path = os.path.join("test", "coassemble", "sra_qc", "SRR3309137_2.fastq.gz")
            self.assertTrue(os.path.exists(sra_r0_path))
            with gzip.open(sra_r0_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SRR3309137.2 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" in file)
                self.assertTrue("@SRR3309137.1 HISEQ06:195:D1DRHACXX:5:1101:1597:2236/1" not in file)

            sra_f1_path = os.path.join("test", "coassemble", "sra_qc", "SRR8334323_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f1_path))
            with gzip.open(sra_f1_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SEQ_ID.1" not in file)

            sra_f2_path = os.path.join("test", "coassemble", "sra_qc", "SRR8334324_1.fastq.gz")
            self.assertTrue(os.path.exists(sra_f2_path))
            with gzip.open(sra_f2_path) as f:
                file = f.readline().decode()
                self.assertTrue("@SEQ_ID.1" in file)

            recover_path = os.path.join("test", "coassemble", "commands", "recover_commands.sh")
            self.assertTrue(os.path.exists(recover_path))
            with open(recover_path) as f:
                file = f.readline()
                self.assertTrue("sra_qc/SRR3309137_1.fastq.gz" in file)
                self.assertTrue("sra_qc/SRR3309137_2.fastq.gz" in file)
                self.assertTrue("sra_qc/SRR8334323_1.fastq.gz" in file)
                self.assertTrue("sra_qc/SRR8334323_2.fastq.gz" in file)
                self.assertTrue("sra_qc/SRR8334324_1.fastq.gz" in file)
                self.assertTrue("sra_qc/SRR8334324_2.fastq.gz" in file)

    def test_update_sra_download_mock_filter_single(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--forward SRR3309137_mismatched SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED_SRA} "
                f"--coassemble-binned {MOCK_BINNED_SRA} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS_SRA_MOCK2} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" --config mock_sra=True\" "
            )
            extern.run(cmd)

            single_ended_path = os.path.join("test", "coassemble", "sra", "single_ended.tsv")
            self.assertTrue(os.path.exists(single_ended_path))
            expected = "\n".join(
                [
                    "\t".join(["sra", "reason"]),
                    "\t".join(["SRR3309137_mismatched", "Consecutive reads do not match (1/10)"]),
                    ""
                ]
            )
            with open(single_ended_path) as f:
                self.assertEqual(expected, f.read())

            elusive_clusters_path = os.path.join("test", "coassemble", "target", "elusive_clusters.tsv")
            self.assertTrue(os.path.exists(elusive_clusters_path))
            expected = "\n".join(
                [
                    "\t".join(["samples", "length", "total_targets", "total_size", "recover_samples", "coassembly"]),
                    "\t".join(["SRR8334324,SRR8334323", "2", "2", "0", "SRR8334323,SRR8334324", "coassembly_0"]),
                    ""
                ]
            )
            with open(elusive_clusters_path) as f:
                self.assertEqual(expected, f.read())

    def test_update_sra_download_mock_fail(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--forward SRR3309137_mismatched SRR8334323 SRR8334324 "
                f"--sra "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED_SRA} "
                f"--coassemble-binned {MOCK_BINNED_SRA} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS_SRA_MOCK3} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" --config mock_sra=True\" "
            )
            with self.assertRaises(Exception):
                extern.run(cmd)

    def test_update_aviary_run(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward SRR8334323 SRR8334324 "
                f"--sra "
                f"--run-aviary "
                f"--aviary-gtdbtk-dir gtdb_release "
                f"--aviary-checkm2-dir CheckM2_database "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED_SRA} "
                f"--coassemble-binned {MOCK_BINNED_SRA} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS_SRA} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-args \" --quiet\" "
            )
            output_comb = extern.run(cmd)

            output_sra = output_comb.split("Building DAG of jobs...")[1]
            self.assertTrue("download_sra" in output_sra)
            self.assertTrue("sra_qc" not in output_sra)
            self.assertTrue("compile_sra_qc" not in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output_sra = output_comb.split("Building DAG of jobs...")[2]
            self.assertTrue("download_sra" not in output_sra)
            self.assertTrue("sra_qc" in output_sra)
            self.assertTrue("compile_sra_qc" in output_sra)
            self.assertTrue("aviary_commands" not in output_sra)

            output = output_comb.split("Building DAG of jobs...")[3]
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("genome_transcripts" not in output)
            self.assertTrue("singlem_pipe_genomes" not in output)
            self.assertTrue("singlem_summarise_genomes" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("query_processing" not in output)
            self.assertTrue("single_assembly" not in output)
            self.assertTrue("count_bp_reads" not in output)
            self.assertTrue("target_elusive" not in output)
            self.assertTrue("cluster_graph" not in output)
            self.assertTrue("collect_genomes" in output)
            self.assertTrue("map_reads" in output)
            self.assertTrue("finish_mapping" in output)
            self.assertTrue("aviary_commands" not in output)
            self.assertTrue("aviary_assemble" in output)
            self.assertTrue("aviary_recover" in output)
            self.assertTrue("aviary_combine" in output)

    def test_update_aviary_dryrun(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--run-aviary "
                f"--aviary-gtdbtk-dir gtdb_release "
                f"--aviary-checkm2-dir CheckM2_database "
                f"--coassemble-output {MOCK_COASSEMBLE} "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
                f"--snakemake-args \" --config aviary_dryrun=True\" "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            coassembly_path = os.path.join("test", "coassemble", "coassemble", "coassembly_0")
            contigs_path = os.path.join(coassembly_path, "assemble", "assembly", "final_contigs.fasta")
            self.assertTrue(os.path.exists(contigs_path))

            assembly_config_path = os.path.join(coassembly_path, "assemble", "config.yaml")
            self.assertTrue(os.path.exists(assembly_config_path))
            assembly_config = load_configfile(assembly_config_path)
            forward_reads = SAMPLE_READS_FORWARD.split(" ")
            self.assertTrue(forward_reads[0] in assembly_config["short_reads_1"])
            self.assertTrue(forward_reads[1] in assembly_config["short_reads_1"])
            self.assertFalse(forward_reads[2] in assembly_config["short_reads_1"])
            self.assertTrue(assembly_config["coassemble"])

            recovery_config_path = os.path.join(coassembly_path, "recover", "config.yaml")
            self.assertTrue(os.path.exists(recovery_config_path))
            recovery_config = load_configfile(recovery_config_path)
            self.assertTrue(forward_reads[0] in recovery_config["short_reads_1"])
            self.assertTrue(forward_reads[1] in recovery_config["short_reads_1"])
            self.assertTrue(forward_reads[2] in recovery_config["short_reads_1"])
            self.assertEqual(recovery_config["fasta"][0], os.path.abspath(contigs_path))
            self.assertEqual(recovery_config["checkm2_db_folder"], "CheckM2_database")
            self.assertEqual(recovery_config["gtdbtk_folder"], "gtdb_release")
            self.assertEqual(recovery_config["eggnog_folder"], ".")

    def test_update_specific_coassembly(self):
        with in_tempdir():
            cmd = (
                f"ibis update "
                f"--assemble-unmapped "
                f"--forward {SAMPLE_READS_FORWARD} "
                f"--reverse {SAMPLE_READS_REVERSE} "
                f"--genomes {GENOMES} "
                f"--coassemble-unbinned {MOCK_UNBINNED} "
                f"--coassemble-binned {MOCK_BINNED} "
                f"--coassemble-targets {MOCK_TARGETS} "
                f"--coassemble-elusive-edges {MOCK_ELUSIVE_EDGES} "
                f"--coassemble-elusive-clusters {MOCK_ELUSIVE_CLUSTERS_TWO} "
                f"--coassemble-summary {MOCK_SUMMARY} "
                f"--coassemblies coassembly_0 "
                f"--output test "
                f"--conda-prefix {path_to_conda} "
            )
            extern.run(cmd)

            config_path = os.path.join("test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))

            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.1.fq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_1_unmapped.2.fq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.1.fq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_2_unmapped.2.fq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_3_unmapped.1.fq.gz")))
            self.assertTrue(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_3_unmapped.2.fq.gz")))
            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_4_unmapped.1.fq.gz")))
            self.assertFalse(os.path.exists(os.path.join("test", "coassemble", "mapping", "sample_4_unmapped.2.fq.gz")))


if __name__ == '__main__':
    unittest.main()
