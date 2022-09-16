#!/usr/bin/env python3

from .__init__ import __version__
__author__ = "Samuel Aroney"

import os
import sys
import logging
import shutil
import subprocess
import bird_tool_utils as btu
from snakemake.io import load_configfile
from ruamel.yaml import YAML

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

    with open(template) as f:
        config = yaml.load(f)

    for key, value in config_items.items():
        if value is not None:
            config[key] = value

    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    logging.info(f"Config file written to: {config_path}")

    return config_path

def copy_input(input, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    logging.info(f"Copying input file {input} to {output}")
    shutil.copyfile(input, output)

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

def cluster(args):
    if not args.forward and not args.sample_singlem:
        raise Exception("Input reads (--forward) or SingleM otu tables (--sample-singlem) must be provided")
    if not args.forward and not args.sample_read_size:
        raise Exception("Input reads (--forward) or read sizes (--sample-read-size) must be provided")
    if args.reverse and not args.forward:
        raise Exception("Reverse reads cannot be provided without forward reads")
    if not args.genome_transcripts and not args.genome_singlem:
        raise Exception("Genome transcripts (--genome-transcripts) or SingleM otu tables (--genome-singlem) must be provided")
    if not args.singlem_metapackage and (not args.sample_singlem or not args.genome_singlem):
        raise Exception("SingleM metapackage (--singlem-metapackage) must be provided when SingleM otu tables not provided")

    output = os.path.abspath(args.output)
    if not os.path.exists(output):
        os.makedirs(output)

    if args.forward:
        forward_reads, reverse_reads = build_reads_list(args.forward, args.reverse)
    if args.sample_read_size:
        with open(args.sample_read_size) as f:
            forward_reads = {line.split(",")[0]: "" for line in f}
            reverse_reads = forward_reads
        copy_input(
            os.path.abspath(args.sample_read_size),
            os.path.join(output, "cluster", os.path.basename(args.sample_read_size)),
        )
    if args.sample_singlem:
        for table in args.sample_singlem:
            copy_input(
                os.path.abspath(table),
                os.path.join(output, "cluster", "summarise", os.path.basename(table))
            )
    if args.genome_transcripts:
        genome_transcripts = {
            os.path.splitext(os.path.basename(transcript))[0]: os.path.abspath(transcript) for transcript in args.genome_transcripts
            }
    if args.genome_singlem:
        copy_input(
            os.path.abspath(args.genome_singlem),
            os.path.join(output, "cluster", "summarise", os.path.basename(args.genome_singlem))
        )

    config_items = {
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "singlem_metapackage": os.path.abspath(args.singlem_metapackage) if args.singlem_metapackage else None,
        "bin_transcripts": genome_transcripts if args.genome_transcripts else None,
        "taxa_of_interest": args.taxa_of_interest if args.taxa_of_interest else None,
        "appraise_sequence_identity": args.appraise_sequence_identity / 100 if args.appraise_sequence_identity > 1 else args.appraise_sequence_identity,
        "min_coassembly_coverage": args.min_sequence_coverage,
        "num_coassembly_samples": args.num_coassembly_samples,
        "max_coassembly_samples": args.max_coassembly_samples if args.max_coassembly_samples else args.num_coassembly_samples,
        "max_coassembly_size": args.max_coassembly_size,
        "max_recovery_samples": args.max_recovery_samples,
    }

    config_path = make_config(
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_cluster.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "cluster.smk",
        output_dir = output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args,
    )

def coassemble(args):
    if not args.forward:
        raise Exception("Input reads must be provided")
    if args.reverse and not args.forward:
        raise Exception("Reverse reads must be provided with forward reads")
    if args.assemble_unmapped and not args.genomes:
        raise Exception("Genomes must be provided for mapping reference with --assemble-unmapped")

    forward_reads, reverse_reads = build_reads_list(args.forward, args.reverse)
    cluster_target_dir = os.path.abspath(os.path.join(args.cluster_output, "target"))

    cluster_appraise_dir = os.path.abspath(os.path.join(args.cluster_output, "appraise"))
    appraise_binned = {}
    for read in forward_reads:
        appraise_binned[read] = os.path.join(cluster_appraise_dir, read + "_binned.otu_table.tsv")
    if args.genomes:
        genomes = {
            os.path.splitext(os.path.basename(genome))[0]: os.path.abspath(genome) for genome in args.genomes
            }

    config_items = {
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "elusive_clusters": os.path.join(cluster_target_dir, "elusive_clusters.tsv"),
        "assemble_unmapped": args.assemble_unmapped,
        "appraise_binned": appraise_binned,
        "genomes": genomes if args.genomes else None,
    }

    output = os.path.abspath(args.output)
    if not os.path.exists(output):
        os.makedirs(output)

    config_path = make_config(
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_coassemble.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "coassemble.smk",
        output_dir = output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args,
    )

def evaluate(args):
    if not args.singlem_metapackage:
        raise Exception("SingleM metapackage must be provided")

    cluster_target_dir = os.path.abspath(os.path.join(args.cluster_output, "target"))
    recovered_bins = {
        os.path.basename(coassembly.rstrip("/")): 
            os.path.abspath(os.path.join(coassembly, "recover", "bins", "final_bins"))
            for coassembly in args.aviary_outputs
        }
    checkm_out = {
        os.path.basename(coassembly.rstrip("/")): 
            os.path.abspath(os.path.join(coassembly, "recover", "bins", "checkm_minimal.tsv"))
            for coassembly in args.aviary_outputs
        }

    config_items = {
        "targets": os.path.join(cluster_target_dir, "targets.tsv"),
        "elusive_edges": os.path.join(cluster_target_dir, "elusive_edges.tsv"),
        "elusive_clusters": os.path.join(cluster_target_dir, "elusive_clusters.tsv"),
        "singlem_metapackage": os.path.abspath(args.singlem_metapackage),
        "recovered_bins": recovered_bins,
        "checkm_out": checkm_out,
        "checkm_version": args.checkm_version,
        "min_completeness": args.min_completeness,
        "max_contamination": args.max_contamination,
    }

    output = os.path.abspath(args.output)
    if not os.path.exists(output):
        os.makedirs(output)

    config_path = make_config(
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_evaluate.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "evaluate.smk",
        output_dir = output,
        cores = args.cores,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
        snakemake_args = args.snakemake_args,
    )

def main():
    main_parser = btu.BirdArgparser(program="Cockatoo", version = __version__,
        examples = {
            "cluster": [
                btu.Example(
                    "cluster reads into suggested coassemblies",
                    "cockatoo cluster --forward reads.1.fq --reverse reads.2.fq --output output_dir"
                ),
                btu.Example(
                    "cluster reads into suggested coassemblies based only on sequences from a specific taxa",
                    "cockatoo cluster --forward reads.1.fq --reverse reads.2.fq --output output_dir --taxa-of-interest \"p__Planctomycetota\""
                ),
                btu.Example(
                    "cluster SingleM outputs (archive otu tables) into suggested coassemblies",
                    "cockatoo cluster --singlem-gzip-archives reads.singlem.json.gz --output output_dir"
                )
            ],
            "coassemble": [
                btu.Example(
                    "coassemble a clustered set of reads",
                    "cockatoo coassemble --cluster-output cluster_dir --output output_dir"
                ),
                btu.Example(
                    "coassemble unmapped reads from a clustered set of reads",
                    "cockatoo coassemble --cluster-output cluster_dir --forward reads.1.fq --reverse reads.2.fq --assemble-unmapped --genomes genome.fna --output output_dir"
                )
            ],
            "evaluate": [
                btu.Example(
                    "evaluate a completed coassembly",
                    "cockatoo evaluate --cluster-output cluster_dir --coassemble-output coassembly_dir --output output_dir"
                )
            ]
        }
        )

    ###########################################################################

    cluster_parser = main_parser.new_subparser("cluster", "Cluster reads into suggested coassemblies by unbinned single-copy marker genes")
    cluster_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    cluster_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    cluster_parser.add_argument("--sample-singlem", nargs='+', help="Summarised SingleM otu tables for each sample. If provided, sample SingleM is skipped")
    cluster_parser.add_argument("--sample-read-size", help="Comma separated list of sample name and size (bp). If provided, sample read counting is skipped")
    cluster_parser.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database")
    cluster_parser.add_argument("--genome-singlem", help="Combined summarised SingleM otu tables for genome transcripts. If provided, genome SingleM is skipped")
    cluster_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    cluster_parser.add_argument("--output", help="Output directory")
    cluster_parser.add_argument("--taxa-of-interest", help="Only consider sequences from this GTDB taxa (e.g. p__Planctomycetota) [default: all]")
    cluster_parser.add_argument("--appraise-sequence-identity", type=int, help="Minimum sequence identity for SingleM appraise against reference database [default: 89%]", default=0.89)
    cluster_parser.add_argument("--min-sequence-coverage", type=int, help="Minimum combined coverage for sequence inclusion [default: 10]", default=10)
    cluster_parser.add_argument("--num-coassembly-samples", type=int, help="Number of samples per coassembly cluster [default: 2]", default=2)
    cluster_parser.add_argument("--max-coassembly-samples", type=int, help="Upper bound for number of samples per coassembly cluster [default: --num-coassembly-samples]", default=None)
    cluster_parser.add_argument("--max-coassembly-size", type=int, help="Maximum size (Gbp) of coassembly cluster [default: None]", default=None)
    cluster_parser.add_argument("--max-recovery-samples", type=int, help="Upper bound for number of related samples to use for differential abundance binning [default: 20]", default=20)
    cluster_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    cluster_parser.add_argument("--cores", type=int, help="Maximum number of cores to use", default=1)
    cluster_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")
    cluster_parser.add_argument("--snakemake-args", help="Additional commands to be supplied to snakemake in the form of a space-prefixed single string e.g. \" --quiet\"", default="")

    ###########################################################################

    coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble suggested coassemblies using Aviary (optionally remove reads mapping to reference genomes)")
    coassemble_parser.add_argument("--cluster-output", help="Output dir from cluster subcommand", required=True)
    coassemble_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    coassemble_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    coassemble_parser.add_argument("--assemble-unmapped", action="store_true", help="Only assemble reads to do not map to reference genomes")
    coassemble_parser.add_argument("--genomes", nargs='+', help="Reference genomes for read mapping")
    coassemble_parser.add_argument("--output", help="output directory")
    coassemble_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    coassemble_parser.add_argument("--cores", type=int, help="Maximum number of cores to use", default=1)
    coassemble_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")
    coassemble_parser.add_argument("--snakemake-args", help="Additional commands to be supplied to snakemake in the form of a space-prefixed single string e.g. \" --quiet\"", default="")

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")
    evaluate_parser.add_argument("--cluster-output", help="Output dir from cluster subcommand", required=True)
    evaluate_parser.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by coassemble subcommand", required=True)
    evaluate_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    evaluate_parser.add_argument("--output", help="output directory")
    evaluate_parser.add_argument("--checkm-version", type=int, help="CheckM version to use to quality cutoffs [default: 2]", default=2)
    evaluate_parser.add_argument("--min-completeness", type=int, help="Include bins with at least this minimum completeness [default: 70]", default=70)
    evaluate_parser.add_argument("--max-contamination", type=int, help="Include bins with at most this maximum contamination [default: 10]", default=10)
    evaluate_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    evaluate_parser.add_argument("--cores", type=int, help="Maximum number of cores to use", default=1)
    evaluate_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")
    evaluate_parser.add_argument("--snakemake-args", help="Additional commands to be supplied to snakemake in the form of a space-prefixed single string e.g. \" --quiet\"", default="")

    ###########################################################################

    args = main_parser.parse_the_args()
    logging.info(f"Cockatoo v{__version__}")
    logging.info(f"Command: {' '.join(['cockatoo'] + sys.argv[1:])}")

    if args.subparser_name == "cluster":
        cluster(args)
    elif args.subparser_name == "coassemble":
        coassemble(args)
    elif args.subparser_name == "evaluate":
        evaluate(args)

if __name__ == "__main__":
    sys.exit(main())
