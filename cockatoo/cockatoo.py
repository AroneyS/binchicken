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

def correlate(args):
    if not args.forward:
        raise Exception("Input reads must be provided")
    if args.reverse and not args.forward:
        raise Exception("Reverse reads must be provided with forward reads")
    if not args.genome_transcripts:
        raise Exception("Genome transcripts must be provided")
    if not args.singlem_metapackage:
        raise Exception("SingleM metapackage must be provided")

    if args.reverse:
        forward_reads = {}
        reverse_reads = {}
        for forward, reverse in zip(args.forward, args.reverse):
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
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_correlate.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "correlate.smk",
        output_dir = output,
        dryrun = args.dryrun,
        conda_prefix = args.conda_prefix,
    )

def evaluate(args):
    if not args.singlem_metapackage:
        raise Exception("SingleM metapackage must be provided")

    correlate_target_dir = os.path.abspath(os.path.join(args.correlate_output, "target"))
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
        "targets": os.path.join(correlate_target_dir, "targets.tsv"),
        "elusive_edges": os.path.join(correlate_target_dir, "elusive_edges.tsv"),
        "elusive_clusters": os.path.join(correlate_target_dir, "elusive_clusters.tsv"),
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
            "correlate": [
                btu.Example(
                    "correlate reads into suggested coassemblies",
                    "cockatoo correlate --forward reads.1.fq --reverse reads.2.fq --output output_dir"
                ),
                btu.Example(
                    "correlate SingleM outputs (archive otu tables) into suggested coassemblies",
                    "cockatoo correlate --singlem-gzip-archives reads.singlem.json.gz --output output_dir"
                )
            ],
            "evaluate": [
                btu.Example(
                    "evaluate a completed coassembly",
                    "cockatoo evaluate --coassembly-output coassembly_dir --output output_dir"
                )
            ]
        }
        )

    ###########################################################################

    correlate_parser = main_parser.new_subparser("correlate", "Correlate reads into suggested coassemblies by unbinned single-copy marker genes")
    correlate_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    correlate_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    correlate_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    correlate_parser.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database")
    correlate_parser.add_argument("--output", help="output directory")
    correlate_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    correlate_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")
    evaluate_parser.add_argument("--correlate-output", help="Output dir from correlate subcommand", required=True)
    evaluate_parser.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by correlate subcommand", required=True)
    evaluate_parser.add_argument("--checkm-version", type=int, help="CheckM version to use for bin evaluation")
    evaluate_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    evaluate_parser.add_argument("--output", help="output directory")
    evaluate_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    evaluate_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

    ###########################################################################

    args = main_parser.parse_the_args()
    logging.info(f"Cockatoo v{__version__}")
    logging.info(f"Command: {' '.join(['cockatoo'] + sys.argv[1:])}")

    if args.subparser_name == "correlate":
        correlate(args)
    elif args.subparser_name == "evaluate":
        evaluate(args)

if __name__ == "__main__":
    sys.exit(main())
