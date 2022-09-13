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
        config[key] = value
    
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    logging.info(f"Config file written to: {config_path}")

    return config_path

def run_workflow(config, workflow, output_dir, cores=16, dryrun=False, snakemake_args="", conda_frontend="mamba"):
    load_configfile(config)

    cmd = (
        "snakemake --snakefile {snakefile} --configfile '{config}' --directory {output_dir} "
        "{jobs} --rerun-incomplete --nolock "
        "--use-conda {conda_frontend} "
        "{dryrun} "
        "{snakemake_args}"
    ).format(
        snakefile=os.path.join(os.path.dirname(__file__), "workflow", workflow),
        config=config,
        output_dir=output_dir,
        jobs=f"--jobs {cores}" if cores is not None else "",
        conda_frontend=f"--conda-frontend {conda_frontend}" if conda_frontend is not None else "",
        dryrun="--dryrun" if dryrun else "",
        snakemake_args=snakemake_args,
    )

    logging.info(f"Executing: {cmd}")
    extern.run(cmd)

def coassemble(args):
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
        os.path.join(os.path.dirname(os.path.realpath(__file__)),"config","template_coassembly.yaml"),
        output,
        config_items
        )

    run_workflow(
        config = config_path,
        workflow = "coassembly.smk",
        output_dir = output,
        dryrun = args.dryrun,
    )

def evaluate(args):
    pass

def main():
    main_parser = btu.BirdArgparser(program="Cockatoo", version = __version__,
        examples = {
            "coassemble": [
                btu.Example(
                    "coassemble reads and bins",
                    "cockatoo coassemble --forward reads.1.fq --reverse reads.2.fq --output output_dir"
                ),
                btu.Example(
                    "coassemble archive otu tables and bins",
                    "cockatoo coassemble --singlem-gzip-archives reads.singlem.json.gz --output output_dir"
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

    coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble reads into contigs and bin")
    coassemble_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    coassemble_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    coassemble_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    coassemble_parser.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database")
    coassemble_parser.add_argument("--output", help="output directory")
    coassemble_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")

    ###########################################################################

    args = main_parser.parse_the_args()
    logging.info(f"Cockatoo v{__version__}")

    if args.subparser_name == "coassemble":
        coassemble(args)
    elif args.subparser_name == "evaluate":
        evaluate(args)

if __name__ == "__main__":
    sys.exit(main())
