#!/usr/bin/env python3

from .__init__ import __version__
__author__ = "Samuel Aroney"

import os
import sys
import logging
import subprocess
import extern
import importlib.resources
import bird_tool_utils as btu
import polars as pl
import polars.selectors as cs
from snakemake.io import load_configfile
from ruamel.yaml import YAML
import copy

FAST_AVIARY_MODE = "fast"
COMPREHENSIVE_AVIARY_MODE = "comprehensive"

def build_reads_list(forward, reverse):
    if reverse:
        if len(forward) != len(reverse):
            raise Exception("Number of forward and reverse reads must be equal")
        forward_reads = {}
        reverse_reads = {}
        for forward, reverse in zip(forward, reverse):
            joint_name = os.path.commonprefix(
                [os.path.basename(forward), os.path.basename(reverse)]
                )
            joint_name = joint_name.rstrip("_").rstrip(".")

            if joint_name in forward_reads:
                raise Exception(f"Duplicate basename: {joint_name}")

            forward_reads[joint_name] = os.path.abspath(forward)
            reverse_reads[joint_name] = os.path.abspath(reverse)
    else:
        forward_reads = {os.path.basename(read): os.path.abspath(read) for read in forward}
        reverse_reads = None

    return forward_reads, reverse_reads

def make_config(template, output_dir, config_items):
    config_path = os.path.join(output_dir, "config.yaml")

    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False

    with importlib.resources.as_file(template) as f:
        config = yaml.load(f)

    for key, value in config_items.items():
        if value is not None:
            config[key] = value

    config = {"Ibis_version": __version__, **config}

    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    logging.info(f"Config file written to: {config_path}")

    return config_path

def load_config(config_path):
    yaml = YAML()
    yaml.version = (1, 1)
    yaml.default_flow_style = False

    with open(config_path) as f:
        config = yaml.load(f)

    return config

def copy_input(input, output, suppress=False):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    if not suppress:
        logging.debug(f"Symbolic-linking input file {input} to {output}")

    try:
        os.symlink(input, output)
    except FileExistsError:
        os.remove(output)
        os.symlink(input, output)

def read_list(path):
    with open(path) as f:
        return [line.strip() for line in f]

def run_workflow(config, workflow, output_dir, cores=16, dryrun=False,
                 snakemake_args="", conda_frontend="mamba", conda_prefix=None):
    load_configfile(config)

    cmd = (
        "snakemake --snakefile {snakefile} --configfile '{config}' --directory {output_dir} "
        "{jobs} --rerun-incomplete --nolock "
        "--use-conda {conda_frontend} {conda_prefix} "
        "{dryrun} "
        "{snakemake_args} "
    ).format(
        snakefile=os.path.join(os.path.dirname(__file__), "workflow", workflow),
        config=config,
        output_dir=output_dir,
        jobs=f"--jobs {cores}" if cores is not None else "",
        conda_frontend=f"--conda-frontend {conda_frontend}" if conda_frontend is not None else "",
        conda_prefix=f"--conda-prefix {conda_prefix}" if conda_prefix is not None else "",
        dryrun="--dryrun" if dryrun else "",
        snakemake_args=snakemake_args,
    )

    logging.info(f"Executing: {cmd}")
    subprocess.check_call(cmd, shell=True)

def download_sra(args):
    config_items = {
        "sra": args.forward,
        "reads_1": {},
        "reads_2": {},
    }

    config_path = make_config(
        importlib.resources.files("ibis.config").joinpath("template_coassemble.yaml"),
        args.output,
        config_items
        )

    if "mock_sra=True" in args.snakemake_args:
        target_rule = "mock_download_sra"
    else:
        target_rule = "download_sra"

    run_workflow(
        config = config_path,
        workflow = "coassemble.smk",
        output_dir = args.output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args + " -- " + target_rule if args.snakemake_args else "-- " + target_rule,
    )

    sra_dir = args.output + "/coassemble/sra/"
    SRA_SUFFIX = ".fastq.gz"
    expected_forward = [sra_dir + f + "_1" + SRA_SUFFIX for f in args.forward]
    expected_reverse = [sra_dir + f + "_2" + SRA_SUFFIX for f in args.forward]

    # Fix single file outputs if interleaved or error if unpaired
    single_ended = {}
    if not args.dryrun and args.sra != "build":
        for f, r in zip(expected_forward, expected_reverse):
            if os.path.isfile(f) & os.path.isfile(r):
                continue

            if os.path.isfile(f) | os.path.isfile(r):
                raise Exception(f"Mismatched downloads with {f} and {r}")

            u = f.replace("_1", "")
            if not os.path.isfile(u):
                raise Exception(f"Missing downloads: {f} and {r}")

            logging.info(f"Checking single-file download for interleaved: {u}")
            cmd = (
                "python {script} "
                "--input {input} "
                "--start-check-pairs {start_check_pairs} "
                "--end-check-pairs {end_check_pairs} "
            ).format(
                script = importlib.resources.files("ibis.workflow.scripts").joinpath("is_interleaved.py"),
                input = u,
                start_check_pairs = 5,
                end_check_pairs = 5,
            )
            output = extern.run(cmd).strip().split("\t")

            if output[0] != "True":
                sra_name = os.path.basename(u).replace(SRA_SUFFIX, "")
                single_ended[sra_name] = output[1]
                logging.warning(f"Download {u} was not interleaved: {output[1]}")
                continue

            logging.info(f"Download {u} was interleaved: {output[1]}")
            logging.info(f"Deinterleaving {u}")
            cmd = (
                "gunzip -c {input} | "
                "paste - - - - - - - - | "
                "tee >(cut -f 1-4 | tr '\t' '\n' | pigz --best --processes {threads} > {output_f}) "
                "| cut -f 5-8 | tr '\t' '\n' | pigz --best --processes {threads} > {output_r}"
            ).format(
                input = u,
                output_f = f,
                output_r = r,
                threads = max(1, args.cores // 2),
            )
            extern.run(cmd)

    if single_ended:
        with open(sra_dir + "single_ended.tsv", "w") as f:
            f.write("sra\treason\n")
            for sra, reason in single_ended.items():
                f.write(f"{sra}\t{reason}\n")

        try:
            elusive_clusters_path = os.path.join(args.output, "coassemble", "target", "elusive_clusters.tsv")
            elusive_clusters = (
                pl.read_csv(elusive_clusters_path, separator="\t")
                .with_columns(
                    single_ended = pl.lit([list(single_ended)]),
                    single_ended_samples =
                        pl.col("samples")
                        .str.split(",")
                        .list.eval(pl.element().is_in(single_ended.keys()))
                        .list.any(),
                    single_ended_recover_samples =
                        pl.col("recover_samples")
                        .str.split(",")
                        .list.eval(pl.element().is_in(single_ended.keys()))
                        .list.any(),
                    )
                .with_columns(
                    pl.col("recover_samples").str.split(",").list.set_difference("single_ended").list.join(","),
                )
            )

            single_ended_sample_coassemblies = (
                elusive_clusters
                .filter(pl.col("single_ended_samples"))
                .get_column("coassembly")
                .to_list()
            )
            for coassembly in single_ended_sample_coassemblies:
                logging.warn(f"Single-ended reads detected in assembly samples for coassembly: {coassembly}, skipping.")

            single_ended_recover_coassemblies = (
                elusive_clusters
                .filter(~pl.col("single_ended_samples"))
                .filter(pl.col("single_ended_recover_samples"))
                .get_column("coassembly")
                .to_list()
            )
            for coassembly in single_ended_recover_coassemblies:
                logging.warn(f"Single-ended reads detected in recovery samples for coassembly: {coassembly}. Removing those samples from recovery.")

            if elusive_clusters.filter(~pl.col("single_ended_samples")).height == 0:
                raise Exception("Single-ended reads detected. All coassemblies contain these reads.")

            os.remove(elusive_clusters_path)
            (
                elusive_clusters
                .filter(~pl.col("single_ended_samples"))
                .select(~cs.starts_with("single_ended"))
                .write_csv(elusive_clusters_path, separator="\t")
            )
        except FileNotFoundError:
            raise Exception("Single-ended reads detected. Remove from input SRA and try again.")

    config_items = {
        "sra": [s for s in args.forward if s not in single_ended],
        "reads_1": {},
        "reads_2": {},
    }

    config_path = make_config(
        importlib.resources.files("ibis.config").joinpath("template_coassemble.yaml"),
        args.output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "coassemble.smk",
        output_dir = args.output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args + " -- compile_sra_qc" if args.snakemake_args else "-- compile_sra_qc",
    )

    # Need to convince Snakemake that the SRA data predates the completed outputs
    # Otherwise it will rerun all the rules with reason "updated input files"
    # Using Jan 1 2000 as pre-date
    sra_qc_dir = args.output + "/coassemble/sra_qc/"
    os.makedirs(sra_qc_dir, exist_ok=True)
    for sra in [f + s for f in args.forward for s in ["_1", "_2"]]:
        subprocess.check_call(f"touch -t 200001011200 {sra_qc_dir + sra + SRA_SUFFIX}", shell=True)

    forward = sorted([sra_qc_dir + f for f in os.listdir(sra_qc_dir) if f.endswith("_1" + SRA_SUFFIX)])
    reverse = sorted([sra_qc_dir + f for f in os.listdir(sra_qc_dir) if f.endswith("_2" + SRA_SUFFIX)])

    return forward, reverse

def set_standard_args(args):
    args.singlem_metapackage = None
    args.sample_singlem = None
    args.sample_singlem_list = None
    args.sample_singlem_dir = None
    args.sample_query = None
    args.sample_query_list = None
    args.sample_query_dir = None
    args.sample_read_size = None
    args.genome_transcripts = None
    args.genome_transcripts_list = None
    args.genome_singlem = None
    args.taxa_of_interest = None
    args.appraise_sequence_identity = 1
    args.min_sequence_coverage = 1
    args.single_assembly = False
    args.no_genomes = False
    args.new_genomes = False
    args.exclude_coassemblies = None
    args.exclude_coassemblies_list = None
    args.num_coassembly_samples = 1
    args.max_coassembly_samples = None
    args.max_coassembly_size = None
    args.max_recovery_samples = 1
    args.prodigal_meta = False

    return(args)

def evaluate_bins(aviary_outputs, checkm_version, min_completeness, max_contamination, iteration=None):
    logging.info(f"Evaluating bins using CheckM{str(checkm_version)} with completeness >= {str(min_completeness)} and contamination <= {str(max_contamination)}")
    recovered_bins = {
        os.path.basename(coassembly.rstrip("/")): 
            os.path.abspath(os.path.join(coassembly, "recover", "bins", "final_bins"))
            for coassembly in aviary_outputs
        }
    checkm_out_dict = {
        os.path.basename(coassembly.rstrip("/")): 
            os.path.abspath(os.path.join(coassembly, "recover", "bins", "checkm_minimal.tsv"))
            for coassembly in aviary_outputs
        }

    if checkm_version == 1:
        completeness_col = "Completeness (CheckM1)"
        contamination_col = "Contamination (CheckM1)"
    elif checkm_version == 2:
        completeness_col = "Completeness (CheckM2)"
        contamination_col = "Contamination (CheckM2)"
    elif checkm_version == "build":
        logging.info("Mock bins for Ibis build")
        return {"iteration_0-coassembly_0-0": os.path.join(aviary_outputs[0], "iteration_0-coassembly_0-0.fna")}
    else:
        raise ValueError("Invalid CheckM version")

    coassembly_bins = {}
    for coassembly in checkm_out_dict:
        checkm_out = pl.read_csv(checkm_out_dict[coassembly], separator="\t")
        passed_bins = checkm_out.filter(
            (pl.col(completeness_col) >= min_completeness) & (pl.col(contamination_col) <= max_contamination),
        ).get_column("Bin Id"
        ).to_list()
        coassembly_bins[coassembly] = passed_bins

    if iteration:
        return {"-".join(["iteration_" + iteration, c, str(i)]): os.path.join(recovered_bins[c], b + ".fna") for c in coassembly_bins for i, b in enumerate(coassembly_bins[c])}
    else:
        return {"-".join([c, str(i)]): os.path.join(recovered_bins[c], b + ".fna") for c in coassembly_bins for i, b in enumerate(coassembly_bins[c])}

def coassemble(args):
    logging.info("Loading sample info")
    if args.forward_list:
        args.forward = read_list(args.forward_list)
    if args.reverse_list:
        args.reverse = read_list(args.reverse_list)
    forward_reads, reverse_reads = build_reads_list(args.forward, args.reverse)

    if args.sample_singlem_list:
        args.sample_singlem = read_list(args.sample_singlem_list)
    if args.sample_singlem:
        for table in args.sample_singlem:
            copy_input(
                os.path.abspath(table),
                os.path.join(args.output, "coassemble", "pipe", os.path.basename(table))
            )
    if args.sample_singlem_dir:
        copy_input(
            os.path.abspath(args.sample_singlem_dir),
            os.path.join(args.output, "coassemble", "pipe")
        )
    if args.sample_query_list:
        args.sample_query = read_list(args.sample_query_list)
    if args.sample_query:
        for table in args.sample_query:
            copy_input(
                os.path.abspath(table),
                os.path.join(args.output, "coassemble", "query", os.path.basename(table))
            )
    if args.sample_query_dir:
        copy_input(
            os.path.abspath(args.sample_query_dir),
            os.path.join(args.output, "coassemble", "query")
        )
    if args.sample_read_size:
        copy_input(
            os.path.abspath(args.sample_read_size),
            os.path.join(args.output, "coassemble", "read_size.csv"),
        )

    if args.exclude_coassemblies_list:
        args.exclude_coassemblies = read_list(args.exclude_coassemblies_list)

    logging.info("Loading genome info")
    if args.genomes_list:
        args.genomes = read_list(args.genomes_list)
    if args.genomes:
        genomes = {
            os.path.splitext(os.path.basename(genome))[0]: os.path.abspath(genome) for genome in args.genomes
            }
    if args.genome_transcripts_list:
        args.genome_transcripts = read_list(args.genome_transcripts_list)
    if args.genome_transcripts:
        for tr in args.genome_transcripts:
            copy_input(
                os.path.abspath(tr),
                os.path.join(args.output, "coassemble", "transcripts", os.path.basename(tr)),
                suppress=True,
            )
    if args.genome_singlem:
        copy_input(
            os.path.abspath(args.genome_singlem),
            os.path.join(args.output, "coassemble", "summarise", "bins_summarised.otu_table.tsv")
        )
    if args.singlem_metapackage:
        metapackage = os.path.abspath(args.singlem_metapackage)
    else:
        try:
            metapackage = os.environ['SINGLEM_METAPACKAGE_PATH']
        except KeyError:
            metapackage = None

    try:
        args.new_genomes
    except AttributeError:
        args.new_genomes = None

    # Load other info
    if args.single_assembly:
        args.num_coassembly_samples = 1
        args.max_coassembly_samples = 1

    config_items = {
        # General config
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "genomes": genomes if args.genomes else None,
        "singlem_metapackage": metapackage,
        # Clustering config
        "taxa_of_interest": args.taxa_of_interest if args.taxa_of_interest else None,
        "appraise_sequence_identity": args.appraise_sequence_identity / 100 if args.appraise_sequence_identity > 1 else args.appraise_sequence_identity,
        "min_coassembly_coverage": args.min_sequence_coverage,
        "single_assembly": args.single_assembly,
        "no_genomes": args.no_genomes,
        "new_genomes": args.new_genomes,
        "exclude_coassemblies": args.exclude_coassemblies,
        "num_coassembly_samples": args.num_coassembly_samples,
        "max_coassembly_samples": args.max_coassembly_samples if args.max_coassembly_samples else args.num_coassembly_samples,
        "max_coassembly_size": args.max_coassembly_size,
        "max_recovery_samples": args.max_recovery_samples,
        "prodigal_meta": args.prodigal_meta,
        # Coassembly config
        "assemble_unmapped": args.assemble_unmapped,
        "unmapping_min_appraised": args.unmapping_min_appraised,
        "unmapping_max_identity": args.unmapping_max_identity,
        "unmapping_max_alignment": args.unmapping_max_alignment,
        "aviary_speed": args.aviary_speed,
        "run_aviary": args.run_aviary,
        "aviary_gtdbtk": args.aviary_gtdbtk_dir,
        "aviary_checkm2": args.aviary_checkm2_dir,
        "aviary_threads": args.aviary_cores,
        "aviary_memory": args.aviary_memory,
        "conda_prefix": args.conda_prefix,
    }

    config_path = make_config(
        importlib.resources.files("ibis.config").joinpath("template_coassemble.yaml"),
        args.output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "coassemble.smk",
        output_dir = args.output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args,
    )

    logging.info(f"Ibis coassemble complete.")
    logging.info(f"Cluster summary at {os.path.join(args.output, 'coassemble', 'summary.tsv')}")
    logging.info(f"More details at {os.path.join(args.output, 'coassemble', 'target', 'elusive_clusters.tsv')}")
    if not args.run_aviary:
        logging.info(f"Aviary commands for coassembly and recovery in shell scripts at {os.path.join(args.output, 'coassemble', 'commands')}")

def evaluate(args):
    logging.info("Loading Ibis coassemble info")
    if args.coassemble_output:
        coassemble_dir = os.path.abspath(args.coassemble_output)
        coassemble_target_dir = os.path.join(coassemble_dir, "target")
        coassemble_appraise_dir = os.path.join(coassemble_dir, "appraise")

        args.coassemble_targets = os.path.join(coassemble_target_dir, "targets.tsv")
        args.coassemble_binned = os.path.join(coassemble_appraise_dir, "binned.otu_table.tsv")
        args.coassemble_elusive_edges = os.path.join(coassemble_target_dir, "elusive_edges.tsv")
        args.coassemble_elusive_clusters = os.path.join(coassemble_target_dir, "elusive_clusters.tsv")
        args.coassemble_summary = os.path.join(coassemble_dir, "summary.tsv")

    if args.new_genomes_list:
        args.new_genomes = read_list(args.new_genomes_list)

    if args.aviary_outputs:
        bins = evaluate_bins(args.aviary_outputs, args.checkm_version, args.min_completeness, args.max_contamination)
    elif args.new_genomes:
        bins = {"-".join([args.coassembly_run, os.path.splitext(os.path.basename(g))[0]]): os.path.abspath(g) for g in args.new_genomes}
    else:
        raise Exception("Programming error: no bins to evaluate")

    if args.singlem_metapackage:
        metapackage = os.path.abspath(args.singlem_metapackage)
    else:
        metapackage = os.environ['SINGLEM_METAPACKAGE_PATH']

    if args.cluster:
        cluster_ani = args.cluster_ani / 100 if args.cluster_ani > 1 else args.cluster_ani
        if args.genomes_list:
            args.genomes = read_list(args.genomes_list)
        original_bins = [os.path.abspath(genome) for genome in args.genomes]
    else:
        cluster_ani = None
        original_bins = None

    config_items = {
        "targets": args.coassemble_targets,
        "binned": args.coassemble_binned,
        "elusive_edges": args.coassemble_elusive_edges,
        "elusive_clusters": args.coassemble_elusive_clusters,
        "coassemble_summary": args.coassemble_summary,
        "singlem_metapackage": metapackage,
        "recovered_bins": bins,
        "checkm_version": args.checkm_version,
        "min_completeness": args.min_completeness,
        "max_contamination": args.max_contamination,
        "cluster": cluster_ani,
        "original_bins": original_bins,
        "prodigal_meta": args.prodigal_meta,
    }

    config_path = make_config(
        importlib.resources.files("ibis.config").joinpath("template_evaluate.yaml"),
        args.output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "evaluate.smk",
        output_dir = args.output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args,
    )

    logging.info(f"Ibis evaluate complete.")
    logging.info(f"Coassembly evaluation summary at {os.path.join(args.output, 'evaluate', 'evaluate', 'summary_stats.tsv')}")
    logging.info(f"Genome recovery breakdown by phyla at {os.path.join(args.output, 'evaluate', 'evaluate', 'plots', 'combined', 'phylum_recovered.png')}")

def update(args):
    logging.info("Loading Ibis coassemble info")
    if args.coassemble_output:
        coassemble_dir = os.path.abspath(args.coassemble_output)
        coassemble_target_dir = os.path.join(coassemble_dir, "target")
        coassemble_appraise_dir = os.path.join(coassemble_dir, "appraise")

        args.coassemble_unbinned = os.path.join(coassemble_appraise_dir, "unbinned.otu_table.tsv")
        args.coassemble_binned = os.path.join(coassemble_appraise_dir, "binned.otu_table.tsv")
        args.coassemble_targets = os.path.join(coassemble_target_dir, "targets.tsv")
        args.coassemble_elusive_edges = os.path.join(coassemble_target_dir, "elusive_edges.tsv")
        args.coassemble_elusive_clusters = os.path.join(coassemble_target_dir, "elusive_clusters.tsv")
        args.coassemble_summary = os.path.join(coassemble_dir, "summary.tsv")

    if args.coassemble_elusive_clusters and not args.coassemblies:
        copy_input(
            os.path.abspath(args.coassemble_elusive_clusters),
            os.path.join(args.output, "coassemble", "target", "elusive_clusters.tsv")
        )
    elif args.coassemble_elusive_clusters and args.coassemblies:
        elusive_clusters = pl.read_csv(os.path.abspath(args.coassemble_elusive_clusters), separator="\t")
        elusive_clusters = elusive_clusters.filter(pl.col("coassembly").is_in(args.coassemblies))

        os.makedirs(os.path.join(args.output, "coassemble", "target"), exist_ok=True)
        elusive_clusters.write_csv(os.path.join(args.output, "coassemble", "target", "elusive_clusters.tsv"), separator="\t")

    copy_input(
        os.path.abspath(args.coassemble_unbinned),
        os.path.join(args.output, "coassemble", "appraise", "unbinned.otu_table.tsv")
    )
    copy_input(
        os.path.abspath(args.coassemble_binned),
        os.path.join(args.output, "coassemble", "appraise", "binned.otu_table.tsv")
    )
    copy_input(
        os.path.abspath(args.coassemble_elusive_edges),
        os.path.join(args.output, "coassemble", "target", "elusive_edges.tsv")
    )
    copy_input(
        os.path.abspath(args.coassemble_targets),
        os.path.join(args.output, "coassemble", "target", "targets.tsv")
    )
    copy_input(
        os.path.abspath(args.coassemble_summary),
        os.path.join(args.output, "coassemble", "summary.tsv")
    )

    if args.sra:
        args.forward, args.reverse = download_sra(args)

    if args.run_aviary:
        args.snakemake_args = args.snakemake_args + " --rerun-triggers mtime -- aviary_combine" if args.snakemake_args else "--rerun-triggers mtime -- aviary_combine"
    else:
        args.snakemake_args = args.snakemake_args + " --rerun-triggers mtime -- aviary_commands" if args.snakemake_args else "--rerun-triggers mtime -- aviary_commands"

    args = set_standard_args(args)
    coassemble(args)

    logging.info(f"Ibis update complete.")
    if not args.run_aviary:
        logging.info(f"Aviary commands for coassembly and recovery in shell scripts at {os.path.join(args.output, 'coassemble', 'commands')}")

def generate_genome_singlem(orig_args, new_genomes):
    args = copy.deepcopy(orig_args)
    args = set_standard_args(args)
    args.singlem_metapackage = orig_args.singlem_metapackage
    args.prodigal_meta = orig_args.prodigal_meta
    args.genomes = new_genomes + ["mock_genome"]

    args.output = os.path.join(orig_args.output, "mock")
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    copy_input(
        os.path.abspath(orig_args.genome_singlem),
        os.path.join(args.output, "coassemble", "pipe", "mock_genome_bin.otu_table.tsv")
        )

    bins_summarised_path = os.path.join(args.output, "coassemble", "summarise", "bins_summarised.otu_table.tsv")
    if args.snakemake_args:
        args.snakemake_args = args.snakemake_args + f" --rerun-triggers mtime -- {bins_summarised_path}"
    else:
        args.snakemake_args = f"--rerun-triggers mtime -- {bins_summarised_path}"

    coassemble(args)

    return os.path.join(args.output, "coassemble", "summarise", "bins_summarised.otu_table.tsv")

def combine_genome_singlem(genome_singlem, new_genome_singlem, path):
    with open(path, "w") as f:
        with open(genome_singlem) as g:
            f.write(g.read())

        with open(new_genome_singlem) as g:
            g.readline()
            f.write(g.read())

def iterate(args):
    logging.info("Evaluating new bins")
    if args.new_genomes_list:
        args.new_genomes = read_list(args.new_genomes_list)

    if args.aviary_outputs:
        bins = evaluate_bins(args.aviary_outputs, args.checkm_version, args.min_completeness, args.max_contamination, args.iteration)
    elif args.new_genomes:
        bins = {os.path.splitext(os.path.basename(g))[0]: os.path.abspath(g) for g in args.new_genomes}
    else:
        raise Exception("Programming error: no bins to evaluate")

    for bin in bins:
        copy_input(
            os.path.abspath(bins[bin]),
            os.path.join(args.output, "recovered_bins", bin + ".fna")
        )
    with open(os.path.join(args.output, "recovered_bins", "bin_provenance.tsv"), "w") as f:
        f.writelines("\n".join(["\t".join([os.path.abspath(bins[bin]), bin + ".fna"]) for bin in bins]))

    new_genomes = [os.path.join(args.output, "recovered_bins", bin + ".fna") for bin in bins]
    args.new_genomes = {
            os.path.splitext(os.path.basename(genome))[0]: os.path.abspath(genome) for genome in new_genomes
            }

    if args.coassemble_output:
        coassemble_appraise_dir = os.path.join(os.path.abspath(args.coassemble_output), "appraise")
        args.coassemble_unbinned = os.path.join(coassemble_appraise_dir, "unbinned.otu_table.tsv")
        args.coassemble_binned = os.path.join(coassemble_appraise_dir, "binned.otu_table.tsv")

    try:
        args.coassemble_unbinned
    except AttributeError:
        args.coassemble_unbinned = None

    try:
        args.coassemble_binned
    except AttributeError:
        args.coassemble_binned = None

    prior_otu_tables = bool(args.coassemble_binned) & bool(args.coassemble_unbinned)

    if not prior_otu_tables and args.genome_singlem and not args.new_genome_singlem:
        args.genome_singlem = generate_genome_singlem(args, new_genomes)
    elif not prior_otu_tables and args.genome_singlem and args.new_genome_singlem:
        combined_genome_singlem = os.path.join(args.output, "coassemble", "pipe", "genome_singlem.otu_table.tsv")
        os.makedirs(os.path.dirname(combined_genome_singlem), exist_ok=True)
        combine_genome_singlem(args.genome_singlem, args.new_genome_singlem, combined_genome_singlem)
        args.genome_singlem = combined_genome_singlem
    elif args.new_genome_singlem:
        copy_input(
            os.path.abspath(args.new_genome_singlem),
            os.path.join(args.output, "coassemble", "summarise", "new_bins_summarised.otu_table.tsv")
        )

    if args.coassemble_unbinned and args.coassemble_binned:
        copy_input(
            os.path.abspath(args.coassemble_unbinned),
            os.path.join(args.output, "coassemble", "appraise", "unbinned_prior.otu_table.tsv")
        )
        copy_input(
            os.path.abspath(args.coassemble_binned),
            os.path.join(args.output, "coassemble", "appraise", "binned_prior.otu_table.tsv")
        )

    if args.genomes_list:
        args.genomes = read_list(args.genomes_list)
        args.genomes_list = None
    args.genomes += new_genomes

    coassemble(args)

    if args.elusive_clusters and not args.dryrun:
        new_cluster = pl.read_csv(os.path.join(args.output, "coassemble", "target", "elusive_clusters.tsv"), separator="\t")
        for cluster in args.elusive_clusters:
            old_cluster = pl.read_csv(cluster, separator="\t")
            comb_cluster = new_cluster.join(old_cluster, on="samples", how="inner")

            if comb_cluster.height > 0:
                _ = comb_cluster.select(
                    pl.col("coassembly").apply(lambda x: logging.warn(f"{x} has been previously suggested"))
                    )
    elif not args.exclude_coassemblies:
        logging.warn("Suggested coassemblies may match those from previous iterations. To check, use `--elusive-clusters`. To exclude manually, use `--exclude-coassembles`")

    logging.info(f"Ibis iterate complete.")
    logging.info(f"Cluster summary at {os.path.join(args.output, 'coassemble', 'summary.tsv')}")
    logging.info(f"More details at {os.path.join(args.output, 'coassemble', 'target', 'elusive_clusters.tsv')}")
    logging.info(f"Aviary commands for coassembly and recovery in shell scripts at {os.path.join(args.output, 'coassemble', 'commands')}")

def build(args):
    output_dir = os.path.join(args.output, "build")
    os.makedirs(output_dir, exist_ok=True)
    conda_prefix = args.conda_prefix
    args.build = True

    # Setup env variables
    os.environ["SNAKEMAKE_CONDA_PREFIX"] = conda_prefix
    extern.run(f"conda env config vars set SNAKEMAKE_CONDA_PREFIX={conda_prefix}")

    # Set args
    args = set_standard_args(args)
    coassemble_config = load_config(importlib.resources.files("ibis.config").joinpath("template_coassemble.yaml"))
    vars(args).update(coassemble_config)
    evaluate_config = load_config(importlib.resources.files("ibis.config").joinpath("template_evaluate.yaml"))
    vars(args).update(evaluate_config)

    args.snakemake_args = args.snakemake_args + " --conda-create-envs-only" if args.snakemake_args else "--conda-create-envs-only"
    args.conda_prefix = conda_prefix

    args.forward_list = None
    args.reverse_list = None
    args.genomes_list = None
    args.new_genomes_list = None
    args.aviary_gtdbtk_dir = "."
    args.aviary_checkm2_dir = "."
    args.aviary_cores = None
    args.assemble_unmapped = True
    args.coassemblies = None
    args.singlem_metapackage = "."

    # Create mock input files
    args.forward = [os.path.join(args.output, "sample_" + s + ".1.fq") for s in ["1", "2", "3"]]
    args.reverse = [os.path.join(args.output, "sample_" + s + ".2.fq") for s in ["1", "2", "3"]]
    args.genomes = [os.path.join(args.output, "genome_1.fna")]

    # Create mock iterate files
    args.checkm_version = "build"
    args.aviary_outputs = [os.path.join(output_dir, "aviary_output")]
    os.makedirs(args.aviary_outputs[0], exist_ok=True)
    new_bins = [os.path.join(args.aviary_outputs[0], "iteration_0-coassembly_0-0.fna")]

    # Create mock update files
    args.coassemble_output = os.path.join(output_dir, "coassemble_output")
    appraise_dir = os.path.join(args.coassemble_output, "appraise")
    os.makedirs(appraise_dir, exist_ok=True)
    otu_tables = [os.path.join(appraise_dir, t + ".otu_table.tsv") for t in ["binned", "unbinned"]]
    target_dir = os.path.join(args.coassemble_output, "target")
    os.makedirs(target_dir, exist_ok=True)
    clusters = [
        os.path.join(target_dir, "elusive_clusters.tsv"),
        os.path.join(target_dir, "elusive_edges.tsv"),
        os.path.join(target_dir, "targets.tsv"),
        ]
    clusters_text = "samples\tlength\ttotal_targets\ttotal_size\trecover_samples\tcoassembly\nSRR8334324,SRR8334323\t2\t2\t0\tSRR8334324,SRR8334323\tcoassembly_0\n"
    with open(clusters[0], "w") as f:
        f.write(clusters_text)

    for item in args.forward + args.reverse + args.genomes + new_bins + otu_tables + clusters:
        subprocess.check_call(f"touch {item}", shell=True)

    logging.info("Building SingleM, CoverM and Prodigal conda environments")
    args.output = os.path.join(output_dir, "build_coassemble")
    os.mkdir(args.output)
    coassemble(args)

    logging.info("Building R conda environments")
    args.output = os.path.join(output_dir, "build_evaluate")
    os.mkdir(args.output)
    evaluate(args)

    logging.info("Building Aviary and Kingfisher conda environments")
    args.output = os.path.join(output_dir, "build_sra")
    os.mkdir(args.output)
    args.run_aviary = True
    args.sra = "build"
    args.forward = ["SRR8334323", "SRR8334324"]
    args.reverse = args.forward
    update(args)

    logging.info(f"Ibis build complete.")
    logging.info(f"Conda envs at {conda_prefix}")

def main():
    main_parser = btu.BirdArgparser(program="Ibis (bin chicken)", version = __version__, program_invocation="ibis",
        examples = {
            "coassemble": [
                btu.Example(
                    "cluster reads into proposed coassemblies",
                    "ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --no-genomes"
                ),
                btu.Example(
                    "cluster reads into proposed coassemblies based on unbinned sequences",
                    "ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
                btu.Example(
                    "cluster reads into proposed coassemblies based on unbinned sequences and coassemble only unbinned reads",
                    "ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --assemble-unmapped"
                ),
                btu.Example(
                    "cluster reads into proposed coassemblies based on unbinned sequences from a specific taxa",
                    "ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --taxa-of-interest \"p__Planctomycetota\""
                ),
                btu.Example(
                    "find relevant samples for differential coverage binning (no coassembly)",
                    "ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --single-assembly"
                ),
            ],
            "evaluate": [
                btu.Example(
                    "evaluate a completed coassembly",
                    "ibis evaluate --coassemble-output coassemble_dir --aviary-outputs coassembly_0_dir ..."
                ),
                btu.Example(
                    "evaluate a completed coassembly by providing genomes directly",
                    "ibis evaluate --coassemble-output coassemble_dir --new-genomes genome_1.fna ... --coassembly-run coassembly_0"
                ),
            ],
            "update": [
                btu.Example(
                    "update previous run to download SRA reads",
                    "ibis update --sra --coassemble-output coassemble_dir --forward SRA000001 ... --genomes genome_1.fna ..."
                ),
                btu.Example(
                    "update previous run to perform unmapping",
                    "ibis update --assemble-unmapped --coassemble-output coassemble_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
                btu.Example(
                    "update previous run to run specific coassemblies",
                    "ibis update --run-aviary --coassemblies coassembly_0 ... --coassemble-output coassemble_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
            ],
            "iterate": [
                btu.Example(
                    "rerun coassemble, adding new bins to database",
                    "ibis iterate --aviary-outputs coassembly_0_dir ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
                btu.Example(
                    "rerun coassemble, adding new bins to database, providing genomes directly",
                    "ibis iterate --new-genomes new_genome_1.fna ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
                btu.Example(
                    "rerun coassemble, adding new bins to database, excluding previous coassembly combinations",
                    "ibis iterate --exclude-coassemblies reads_1,reads_2 --aviary-outputs coassembly_0_dir ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ..."
                ),
            ],
            "build": [
                btu.Example(
                    "create dependency conda environments",
                    "ibis build --conda-prefix path_to_conda"
                ),
            ],
        }
        )

    ###########################################################################
    def add_general_snakemake_options(argument_group, required_conda_prefix=False):
        argument_group.add_argument("--output", help="Output directory [default: .]", default="./")
        argument_group.add_argument("--conda-prefix", help="Path to conda environment install location. default: Use path from CONDA_ENV_PATH env variable", default=None, required=required_conda_prefix)
        argument_group.add_argument("--cores", type=int, help="Maximum number of cores to use", default=1)
        argument_group.add_argument("--dryrun", action="store_true", help="dry run workflow")
        argument_group.add_argument("--snakemake-args", help="Additional commands to be supplied to snakemake in the form of a space-prefixed single string e.g. \" --quiet\"", default="")
    
    def add_base_arguments(argument_group):
        argument_group.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
        argument_group.add_argument("--forward-list", "--reads-list", "--sequences-list", help="input forward/unpaired nucleotide read sequence(s) newline separated")
        argument_group.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
        argument_group.add_argument("--reverse-list", help="input reverse nucleotide read sequence(s) newline separated")
        argument_group.add_argument("--genomes", nargs='+', help="Reference genomes for read mapping")
        argument_group.add_argument("--genomes-list", help="Reference genomes for read mapping newline separated")

    def add_evaluation_options(argument_group):
        argument_group.add_argument("--checkm-version", type=int, help="CheckM version to use to quality cutoffs [default: 2]", default=2)
        argument_group.add_argument("--min-completeness", type=int, help="Include bins with at least this minimum completeness [default: 70]", default=70)
        argument_group.add_argument("--max-contamination", type=int, help="Include bins with at most this maximum contamination [default: 10]", default=10)

    def add_aviary_options(argument_group):
        argument_group.add_argument("--aviary-speed", help="Run Aviary recover in 'fast' or 'comprehensive' mode. Fast mode skips slow binners and refinement steps.",
                                    default=FAST_AVIARY_MODE, choices=[FAST_AVIARY_MODE, COMPREHENSIVE_AVIARY_MODE])
        argument_group.add_argument("--run-aviary", action="store_true", help="Run Aviary commands for all identified coassemblies (unless specified)")
        argument_group.add_argument("--aviary-gtdbtk-dir", help="Path to GTDB-Tk database directory for Aviary")
        argument_group.add_argument("--aviary-checkm2-dir", help="Path to CheckM2 database directory for Aviary")
        argument_group.add_argument("--aviary-cores", type=int, help="Maximum number of cores for Aviary to use. Half used for recovery.", default=64)
        argument_group.add_argument("--aviary-memory", type=int, help="Maximum amount of memory for Aviary to use (Gigabytes). Half used for recovery", default=500)

    def add_main_coassemble_output_arguments(argument_group):
        argument_group.add_argument("--coassemble-output", help="Output dir from cluster subcommand")
        argument_group.add_argument("--coassemble-unbinned", help="SingleM appraise unbinned output from Ibis coassemble (alternative to --coassemble-output)")
        argument_group.add_argument("--coassemble-binned", help="SingleM appraise binned output from Ibis coassemble (alternative to --coassemble-output)")

    def add_coassemble_output_arguments(argument_group):
        add_main_coassemble_output_arguments(argument_group)
        argument_group.add_argument("--coassemble-targets", help="Target sequences output from Ibis coassemble (alternative to --coassemble-output)")
        argument_group.add_argument("--coassemble-elusive-edges", help="Elusive edges output from Ibis coassemble (alternative to --coassemble-output)")
        argument_group.add_argument("--coassemble-elusive-clusters", help="Elusive clusters output from Ibis coassemble (alternative to --coassemble-output)")
        argument_group.add_argument("--coassemble-summary", help="Summary output from Ibis coassemble (alternative to --coassemble-output)")

    ###########################################################################

    coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble reads clustered by unbinned single-copy marker genes")

    def add_coassemble_arguments(parser):
        # Base arguments
        coassemble_base = parser.add_argument_group("Base input arguments")
        add_base_arguments(coassemble_base)
        coassemble_base.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
        # Midpoint arguments
        coassemble_midpoint = parser.add_argument_group("Intermediate results input arguments")
        coassemble_midpoint.add_argument("--sample-singlem", nargs='+', help="SingleM otu tables for each sample, in the form \"[sample name]_read.otu_table.tsv\". If provided, SingleM pipe sample is skipped")
        coassemble_midpoint.add_argument("--sample-singlem-list", help="SingleM otu tables for each sample, in the form \"[sample name]_read.otu_table.tsv\" newline separated. If provided, SingleM pipe sample is skipped")
        coassemble_midpoint.add_argument("--sample-singlem-dir", help="Directory containing SingleM otu tables for each sample, in the form \"[sample name]_read.otu_table.tsv\". If provided, SingleM pipe sample is skipped")
        coassemble_midpoint.add_argument("--sample-query", nargs='+', help="Queried SingleM otu tables for each sample against genome database, in the form \"[sample name]_query.otu_table.tsv\". If provided, SingleM pipe and appraise are skipped")
        coassemble_midpoint.add_argument("--sample-query-list", help="Queried SingleM otu tables for each sample against genome database, in the form \"[sample name]_query.otu_table.tsv\" newline separated. If provided, SingleM pipe and appraise are skipped")
        coassemble_midpoint.add_argument("--sample-query-dir", help="Directory containing Queried SingleM otu tables for each sample against genome database, in the form \"[sample name]_query.otu_table.tsv\". If provided, SingleM pipe and appraise are skipped")
        coassemble_midpoint.add_argument("--sample-read-size", help="Comma separated list of sample name and size (bp). If provided, sample read counting is skipped")
        coassemble_midpoint.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database, in the form \"[genome]_protein.fna\"")
        coassemble_midpoint.add_argument("--genome-transcripts-list", help="Genome transcripts for reference database, in the form \"[genome]_protein.fna\" newline separated")
        coassemble_midpoint.add_argument("--genome-singlem", help="Combined SingleM otu tables for genome transcripts. If provided, genome SingleM is skipped")
        # Clustering options
        coassemble_clustering = parser.add_argument_group("Clustering options")
        coassemble_clustering.add_argument("--taxa-of-interest", help="Only consider sequences from this GTDB taxa (e.g. p__Planctomycetota) [default: all]")
        coassemble_clustering.add_argument("--appraise-sequence-identity", type=int, help="Minimum sequence identity for SingleM appraise against reference database [default: 86%, Genus-level]", default=0.86)
        coassemble_clustering.add_argument("--min-sequence-coverage", type=int, help="Minimum combined coverage for sequence inclusion [default: 10]", default=10)
        coassemble_clustering.add_argument("--no-genomes", action="store_true", help="Run pipeline without genomes")
        coassemble_clustering.add_argument("--single-assembly", action="store_true", help="Skip appraise to discover samples to differential abundance binning. Forces --num-coassembly-samples and --max-coassembly-samples to 1")
        coassemble_clustering.add_argument("--exclude-coassemblies", nargs='+', help="List of coassemblies to exclude, space separated, in the form \"sample_1,sample_2\"")
        coassemble_clustering.add_argument("--exclude-coassemblies-list", help="List of coassemblies to exclude, space separated, in the form \"sample_1,sample_2\", newline separated")
        coassemble_clustering.add_argument("--num-coassembly-samples", type=int, help="Number of samples per coassembly cluster [default: 2]", default=2)
        coassemble_clustering.add_argument("--max-coassembly-samples", type=int, help="Upper bound for number of samples per coassembly cluster [default: --num-coassembly-samples]", default=None)
        coassemble_clustering.add_argument("--max-coassembly-size", type=int, help="Maximum size (Gbp) of coassembly cluster [default: None]", default=None)
        coassemble_clustering.add_argument("--max-recovery-samples", type=int, help="Upper bound for number of related samples to use for differential abundance binning [default: 20]", default=20)
        coassemble_clustering.add_argument("--prodigal-meta", action="store_true", help="Use prodigal \"-p meta\" argument (for testing)")
        # Coassembly options
        coassemble_coassembly = parser.add_argument_group("Coassembly options")
        coassemble_coassembly.add_argument("--assemble-unmapped", action="store_true", help="Only assemble reads that do not map to reference genomes")
        coassemble_coassembly.add_argument("--unmapping-min-appraised", type=int, help="Minimum fraction of sequences binned to justify unmapping [default: 0.1]", default=0.1)
        coassemble_coassembly.add_argument("--unmapping-max-identity", type=float, help="Maximum sequence identity of mapped sequences kept for coassembly [default: 99%]", default=99)
        coassemble_coassembly.add_argument("--unmapping-max-alignment", type=float, help="Maximum percent alignment of mapped sequences kept for coassembly [default: 99%]", default=99)
        add_aviary_options(coassemble_coassembly)
        # General options
        coassemble_general = parser.add_argument_group("General options")
        add_general_snakemake_options(coassemble_general)

    add_coassemble_arguments(coassemble_parser)

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")
    # Base arguments
    evaluate_base = evaluate_parser.add_argument_group("Base input arguments")
    add_coassemble_output_arguments(evaluate_base)
    evaluate_base.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by coassemble subcommand")
    evaluate_base.add_argument("--new-genomes", nargs='+', help="New genomes to evaluate (alternative to --aviary-outputs, also requires --coassembly-run)")
    evaluate_base.add_argument("--new-genomes-list", help="New genomes to evaluate (alternative to --aviary-outputs, also requires --coassembly-run) newline separated")
    evaluate_base.add_argument("--coassembly-run", help="Name of coassembly run to produce new genomes (alternative to --aviary-outputs, also requires --new-genomes)")
    evaluate_base.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    evaluate_base.add_argument("--prodigal-meta", action="store_true", help="Use prodigal \"-p meta\" argument (for testing)")
    # Evaluate options
    evaluate_evaluation = evaluate_parser.add_argument_group("Evaluation options")
    add_evaluation_options(evaluate_evaluation)
    # Cluster options
    evaluate_cluster = evaluate_parser.add_argument_group("Cluster options")
    evaluate_cluster.add_argument("--cluster", action="store_true", help="Cluster new and original genomes and report number of new clusters")
    evaluate_cluster.add_argument("--cluster-ani", type=int, help="Cluster using this sequence identity [default: 86%]", default=0.86)
    evaluate_cluster.add_argument("--genomes", nargs='+', help="Original genomes used as references for coassemble subcommand")
    evaluate_cluster.add_argument("--genomes-list", help="Original genomes used as references for coassemble subcommand newline separated")
    # General options
    evaluate_general = evaluate_parser.add_argument_group("General options")
    add_general_snakemake_options(evaluate_general)

    ###########################################################################

    update_parser = main_parser.new_subparser("update", "Coassemble pre-clustered reads")
    # Base arguments
    update_base = update_parser.add_argument_group("Input arguments")
    add_base_arguments(update_base)
    update_base.add_argument("--sra", action="store_true", help="Download reads from SRA (read argument still required)")
    # Coassembly options
    update_coassembly = update_parser.add_argument_group("Coassembly options")
    add_coassemble_output_arguments(update_coassembly)
    update_coassembly.add_argument("--coassemblies", nargs='+', help="Choose specific coassemblies from elusive clusters (e.g. coassembly_0)")
    update_coassembly.add_argument("--assemble-unmapped", action="store_true", help="Only assemble reads that do not map to reference genomes")
    update_coassembly.add_argument("--unmapping-min-appraised", type=float, help="Minimum fraction of sequences binned to justify unmapping [default: 0.1]", default=0.1)
    update_coassembly.add_argument("--unmapping-max-identity", type=float, help="Maximum sequence identity of mapped sequences kept for coassembly [default: 99%]", default=99)
    update_coassembly.add_argument("--unmapping-max-alignment", type=float, help="Maximum percent alignment of mapped sequences kept for coassembly [default: 99%]", default=99)
    add_aviary_options(update_coassembly)
    # General options
    update_general = update_parser.add_argument_group("General options")
    add_general_snakemake_options(update_general)

    ###########################################################################

    iterate_parser = main_parser.new_subparser("iterate", "Iterate coassemble using new bins")
    # Iterate options
    iterate_iteration = iterate_parser.add_argument_group("Iteration options")
    iterate_iteration.add_argument("--iteration", help="Iteration number used for unique bin naming", default="0")
    iterate_iteration.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by coassemble subcommand")
    iterate_iteration.add_argument("--new-genomes", nargs='+', help="New genomes to iterate (alternative to --aviary-outputs)")
    iterate_iteration.add_argument("--new-genomes-list", help="New genomes to iterate (alternative to --aviary-outputs) newline separated")
    iterate_iteration.add_argument("--new-genome-singlem", help="Combined SingleM otu tables for new genome transcripts. If provided, genome SingleM is skipped")
    iterate_iteration.add_argument("--elusive-clusters", nargs='+', help="Previous elusive_clusters.tsv files produced by coassemble subcommand (used to check for duplicated coassembly suggestions)")
    add_main_coassemble_output_arguments(iterate_iteration)
    add_evaluation_options(iterate_iteration)
    # Coassembly options
    add_coassemble_arguments(iterate_parser)

    ###########################################################################

    build_parser = main_parser.new_subparser("build", "Create dependency conda environments")
    add_general_snakemake_options(build_parser, required_conda_prefix=True)

    ###########################################################################

    args = main_parser.parse_the_args()
    logging.info(f"Ibis v{__version__}")
    logging.info(f"Command: {' '.join(['ibis'] + sys.argv[1:])}")

    args.output = os.path.abspath(args.output)
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    def base_argument_verification(args):
        if not args.forward and not args.forward_list:
            raise Exception("Input reads must be provided")
        if not args.reverse and not args.reverse_list:
            try:
                if args.sra:
                    logging.info("SRA reads reverse reads not required")
                else:
                    raise Exception("Interleaved and long-reads not yet implemented")
            except AttributeError:
                raise Exception("Interleaved and long-reads not yet implemented")
        if not (args.genomes or args.genomes_list or args.no_genomes):
            raise Exception("Input genomes must be provided, or argument --no-genomes must be used")
        if (args.forward and args.forward_list) or (args.reverse and args.reverse_list) or (args.genomes and args.genomes_list):
            raise Exception("General and list arguments are mutually exclusive")

    def coassemble_argument_verification(args):
        base_argument_verification(args)
        if (args.sample_query or args.sample_query_list or args.sample_query_dir) and not (args.sample_singlem or args.sample_singlem_list or args.sample_singlem_dir):
            raise Exception("Input SingleM query (--sample-query) requires SingleM otu tables (--sample-singlem) for coverage")
        if args.assemble_unmapped and args.single_assembly:
            raise Exception("Assemble unmapped is incompatible with single-sample assembly")
        if args.assemble_unmapped and not args.genomes and not args.genomes_list:
            raise Exception("Reference genomes must be provided to assemble unmapped reads")
        if not args.singlem_metapackage and not os.environ['SINGLEM_METAPACKAGE_PATH'] and not args.sample_query and not args.sample_query_list:
            raise Exception("SingleM metapackage (--singlem-metapackage or SINGLEM_METAPACKAGE_PATH environment variable, see SingleM data) must be provided when SingleM query otu tables are not provided")
        if (args.sample_singlem and args.sample_singlem_list) or (args.sample_singlem_dir and args.sample_singlem_list) or (args.sample_singlem and args.sample_singlem_dir) or \
            (args.sample_query and args.sample_query_list) or (args.sample_query_dir and args.sample_query_list) or (args.sample_query and args.sample_query_dir):
            raise Exception("General, list and directory arguments are mutually exclusive")
        if args.max_coassembly_samples:
            if args.max_coassembly_samples > args.max_recovery_samples:
                raise Exception("Max recovery samples (--max-recovery-samples) must be greater than or equal to max coassembly samples (--max-coassembly-samples)")
        else:
            if args.num_coassembly_samples > args.max_recovery_samples:
                raise Exception("Max recovery samples (--max-recovery-samples) must be greater than or equal to number of coassembly samples (--num-coassembly-samples)")
        if args.run_aviary and not (args.aviary_gtdbtk_dir and args.aviary_checkm2_dir):
            raise Exception("Run Aviary (--run-aviary) requires paths to GTDB-Tk and CheckM2 databases to be provided (--aviary-gtdbtk-dir and --aviary-checkm2-dir)")
        if (args.sample_query or args.sample_query_list or args.sample_query_dir) and args.taxa_of_interest and args.assemble_unmapped:
            raise Exception("Unmapping is incompatible with the combination of sample query and taxa of interest")

    def coassemble_output_argument_verification(args):
        if not args.coassemble_output and not (args.coassemble_unbinned and args.coassemble_binned and args.coassemble_targets and \
                                               args.coassemble_elusive_edges and args.coassemble_elusive_clusters and args.coassemble_summary):
            raise Exception("Either Ibis coassemble output (--coassemble-output) or specific input files must be provided")

    if args.subparser_name == "coassemble":
        coassemble_argument_verification(args)
        coassemble(args)

    elif args.subparser_name == "evaluate":
        coassemble_output_argument_verification(args)
        if not args.singlem_metapackage and not os.environ['SINGLEM_METAPACKAGE_PATH']:
            raise Exception("SingleM metapackage (--singlem-metapackage or SINGLEM_METAPACKAGE_PATH environment variable, see SingleM data) must be provided")
        if args.cluster and not (args.genomes or args.genomes_list):
            raise Exception("Reference genomes must be provided to cluster with new genomes")
        if not args.aviary_outputs and not (args.new_genomes or args.new_genomes_list):
            raise Exception("New genomes or aviary outputs must be provided for evaluation")
        if (args.new_genomes or args.new_genomes_list) and not args.coassembly_run:
            raise Exception("Name of coassembly run must be provided to evaluate binning")
        evaluate(args)

    elif args.subparser_name == "update":
        base_argument_verification(args)
        coassemble_output_argument_verification(args)
        if args.run_aviary and not (args.aviary_gtdbtk_dir and args.aviary_checkm2_dir):
            raise Exception("Run Aviary (--run-aviary) requires paths to GTDB-Tk and CheckM2 databases to be provided (--aviary-gtdbtk-dir and --aviary-checkm2-dir)")
        update(args)

    elif args.subparser_name == "iterate":
        if args.sample_query or args.sample_query_list or args.sample_query_dir:
            raise Exception("Query arguments are incompatible with Ibis iterate")
        if args.sample_singlem_dir or args.sample_query_dir:
            raise Exception("Directory arguments are incompatible with Ibis iterate")
        if args.single_assembly:
            raise Exception("Single assembly is incompatible with Ibis iterate")
        if not args.aviary_outputs and not (args.new_genomes or args.new_genomes_list):
            raise Exception("New genomes or aviary outputs must be provided for iteration")
        coassemble_argument_verification(args)
        iterate(args)

    elif args.subparser_name == "build":
        build(args)

if __name__ == "__main__":
    sys.exit(main())
