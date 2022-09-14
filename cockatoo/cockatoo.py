#!/usr/bin/env python3

from .__init__ import __version__
__author__ = "Samuel Aroney"

import os
import sys
import logging
import bird_tool_utils as btu
import extern
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
            joint_name = joint_name.rstrip(".")

            if joint_name in forward_reads:
                raise Exception(f"Duplicate basename: {joint_name}")

            forward_reads[joint_name] = os.path.abspath(forward)
            reverse_reads[joint_name] = os.path.abspath(reverse)
    else:
        forward_reads = {os.path.basename(read): os.path.abspath(read) for read in args.forward}
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
    extern.run(cmd)

def cluster(args):
    if not args.forward:
        raise Exception("Input reads must be provided")
    if args.reverse and not args.forward:
        raise Exception("Reverse reads must be provided with forward reads")
    if not args.genome_transcripts:
        raise Exception("Genome transcripts must be provided")
    if not args.singlem_metapackage:
        raise Exception("SingleM metapackage must be provided")

    forward_reads, reverse_reads = build_reads_list(args.forward, args.reverse)
    if args.genome_transcripts:
        genome_transcripts = {
            os.path.splitext(os.path.basename(transcript))[0]: os.path.abspath(transcript) for transcript in args.genome_transcripts
            }

    config_items = {
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "singlem_metapackage": os.path.abspath(args.singlem_metapackage),
        "bin_transcripts": genome_transcripts if args.genome_transcripts else None,
    }

    output = os.path.abspath(args.output)
    if not os.path.exists(output):
        os.makedirs(output)

    config_path = make_config(
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_cluster.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "cluster.smk",
        output_dir = output,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
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
    if args.genomes:
        genomes = {
            os.path.splitext(os.path.basename(genome))[0]: os.path.abspath(genome) for genome in args.genomes
            }

    config_items = {
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "elusive_clusters": os.path.join(cluster_target_dir, "elusive_clusters.tsv"),
        "assemble_unmapped": args.assemble_unmapped,
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
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
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
        "checkm_version": args.checkm_version if args.checkm_version else None,
        "recovered_bins": recovered_bins,
        "checkm_out": checkm_out,
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
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
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
                    "cluster SingleM outputs (archive otu tables) into suggested coassemblies",
                    "cockatoo cluster --singlem-gzip-archives reads.singlem.json.gz --output output_dir"
                )
            ],
            "coassemble": [
                btu.Example(
                    "coassemble a clustered set of reads",
                    "cockatoo coassemble --cluster-output coassembly_dir --output output_dir"
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
    cluster_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    cluster_parser.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database")
    cluster_parser.add_argument("--output", help="output directory")
    cluster_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    cluster_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

    ###########################################################################

    coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble suggested coassemblies using Aviary (optionally remove reads mapping to reference genomes)")
    coassemble_parser.add_argument("--cluster-output", help="Output dir from cluster subcommand", required=True)
    coassemble_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    coassemble_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    coassemble_parser.add_argument("--assemble-unmapped", action="store_true", help="Only assemble reads to do not map to reference genomes")
    coassemble_parser.add_argument("--genomes", nargs='+', help="Reference genomes for read mapping")
    coassemble_parser.add_argument("--output", help="output directory")
    coassemble_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    coassemble_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")
    evaluate_parser.add_argument("--cluster-output", help="Output dir from cluster subcommand", required=True)
    evaluate_parser.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by coassemble subcommand", required=True)
    evaluate_parser.add_argument("--checkm-version", type=int, help="CheckM version to use for bin evaluation")
    evaluate_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    evaluate_parser.add_argument("--output", help="output directory")
    evaluate_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    evaluate_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

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
