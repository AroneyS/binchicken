#!/usr/bin/env python3

from .__init__ import __version__
__author__ = "Samuel Aroney"

import os
import sys
import logging
import shutil
import subprocess
import importlib.resources
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

    with importlib.resources.as_file(template) as f:
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

def coassemble(args):
    if not args.forward and not args.forward_list:
        raise Exception("Input reads must be provided")
    if args.sample_query and not args.sample_singlem:
        raise Exception("Input SingleM query (--sample-query) requires SingleM otu tables (--sample-singlem) for coverage")
    if not args.genomes and not args.genomes_list and not args.single_assembly:
        raise Exception("Reference genomes must be provided")
    if not args.singlem_metapackage and not args.sample_query:
        raise Exception("SingleM metapackage (--singlem-metapackage) must be provided when SingleM query otu tables are not provided")
    if (args.forward and args.forward_list) or (args.reverse and args.reverse_list) or (args.genomes and args.genomes_list):
        raise Exception("General argument cannot be provided with list argument")

    output = os.path.abspath(args.output)
    if not os.path.exists(output):
        os.makedirs(output)

    # Load sample info
    if args.forward_list:
        args.forward = read_list(args.forward_list)
    if args.reverse_list:
        args.reverse = read_list(args.reverse_list)
    forward_reads, reverse_reads = build_reads_list(args.forward, args.reverse)
    if args.sample_singlem:
        for table in args.sample_singlem:
            copy_input(
                os.path.abspath(table),
                os.path.join(output, "coassemble", "pipe", os.path.basename(table))
            )
    if args.sample_query:
        for table in args.sample_query:
            copy_input(
                os.path.abspath(table),
                os.path.join(output, "coassemble", "query", os.path.basename(table))
            )
    if args.sample_read_size:
        copy_input(
            os.path.abspath(args.sample_read_size),
            os.path.join(output, "coassemble", "read_size.csv"),
        )

    # Load genome info
    if args.genomes_list:
        args.genomes = read_list(args.genomes_list)
    if args.genomes:
        genomes = {
            os.path.splitext(os.path.basename(genome))[0]: os.path.abspath(genome) for genome in args.genomes
            }
    if args.genome_transcripts_list:
        args.genome_transcripts = read_list(args.genome_transcripts_list)
    if args.genome_transcripts:
        genome_transcripts = {
            os.path.splitext(os.path.basename(transcript))[0]: os.path.abspath(transcript) for transcript in args.genome_transcripts
            }
    if args.genome_singlem:
        copy_input(
            os.path.abspath(args.genome_singlem),
            os.path.join(output, "coassemble", "summarise", "bins_summarised.otu_table.tsv")
        )

    # Load other info
    if args.single_assembly:
        args.num_coassembly_samples = 1
        args.max_coassembly_samples = 1

    config_items = {
        # Sample config
        "reads_1": forward_reads,
        "reads_2": reverse_reads,
        "singlem_metapackage": os.path.abspath(args.singlem_metapackage) if args.singlem_metapackage else None,
        # Genome config
        "genomes": genomes if args.genomes else None,
        "bin_transcripts": genome_transcripts if args.genome_transcripts else None,
        # Clustering config
        "taxa_of_interest": args.taxa_of_interest if args.taxa_of_interest else None,
        "appraise_sequence_identity": args.appraise_sequence_identity / 100 if args.appraise_sequence_identity > 1 else args.appraise_sequence_identity,
        "min_coassembly_coverage": args.min_sequence_coverage,
        "single_assembly": args.single_assembly,
        "num_coassembly_samples": args.num_coassembly_samples,
        "max_coassembly_samples": args.max_coassembly_samples if args.max_coassembly_samples else args.num_coassembly_samples,
        "max_coassembly_size": args.max_coassembly_size,
        "max_recovery_samples": args.max_recovery_samples,
        # Coassembly config
        "assemble_unmapped": args.assemble_unmapped,
        "abstract_options": args.abstract_options,
        "aviary_threads": args.aviary_cores,
        "aviary_memory": args.aviary_memory,
    }

    config_path = make_config(
        importlib.resources.files("cockatoo.config").joinpath("template_coassemble.yaml"),
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
        importlib.resources.files("cockatoo.config").joinpath("template_evaluate.yaml"),
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
                    "cluster reads into suggested coassemblies based on unbinned sequences",
                    "cockatoo cluster --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genome-transcripts genome_protein.fna ... --singlem-metapackage metapackage.smpkg"
                ),
                btu.Example(
                    "cluster reads into suggested coassemblies based on unbinned sequences from a specific taxa",
                    "cockatoo cluster --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genome-transcripts genome_protein.fna ... --taxa-of-interest \"p__Planctomycetota\" --singlem-metapackage metapackage.smpkg"
                ),
                btu.Example(
                    "cluster SingleM outputs into suggested coassemblies (skips SingleM pipe)",
                    "cockatoo cluster --sample-singlem reads_1.otu_table.tsv ... --sample-read-size read_size.csv --genome-singlem genome.otu_table.tsv --singlem-metapackage metapackage.smpkg"
                ),
                btu.Example(
                    "cluster SingleM query outputs into suggested coassemblies (skips SingleM pipe and appraise)",
                    "cockatoo cluster --sample-query reads_1_query.otu_table.tsv ... --sample-read-size read_size.csv"
                ),
                btu.Example(
                    "find relevant samples for differential coverage binning (no coassembly)",
                    "cockatoo cluster --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --single-assembly --singlem-metapackage metapackage.smpkg"
                ),
            ],
            "coassemble": [
                btu.Example(
                    "coassemble a clustered set of reads",
                    "cockatoo coassemble --cluster-output cluster_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ..."
                ),
                btu.Example(
                    "coassemble unmapped reads from a clustered set of reads",
                    "cockatoo coassemble --cluster-output cluster_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome.fna ... --assemble-unmapped"
                ),
            ],
            "evaluate": [
                btu.Example(
                    "evaluate a completed coassembly",
                    "cockatoo evaluate --cluster-output cluster_dir --aviary-outputs coassembly_0_dir ... --singlem-metapackage metapackage.smpkg"
                ),
            ]
        }
        )

    ###########################################################################

    coassemble_parser = main_parser.new_subparser("coassemble", "Coassemble reads clustered by unbinned single-copy marker genes")
    # Base arguments
    coassemble_parser.add_argument("--forward", "--reads", "--sequences", nargs='+', help="input forward/unpaired nucleotide read sequence(s)")
    coassemble_parser.add_argument("--forward-list", "--reads-list", "--sequences-list", help="input forward/unpaired nucleotide read sequence(s) newline separated")
    coassemble_parser.add_argument("--reverse", nargs='+', help="input reverse nucleotide read sequence(s)")
    coassemble_parser.add_argument("--reverse-list", help="input reverse nucleotide read sequence(s) newline separated")
    coassemble_parser.add_argument("--genomes", nargs='+', help="Reference genomes for read mapping")
    coassemble_parser.add_argument("--genomes-list", help="Reference genomes for read mapping newline separated")
    coassemble_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    # Midpoint arguments
    coassemble_parser.add_argument("--sample-singlem", nargs='+', help="SingleM otu tables for each sample, in the form \"[sample name]_read.otu_table.tsv\". If provided, SingleM pipe sample is skipped")
    coassemble_parser.add_argument("--sample-query", nargs='+', help="Queried SingleM otu tables for each sample against genome database, in the form \"[sample name]_query.otu_table.tsv\". If provided, SingleM pipe and appraise are skipped")
    coassemble_parser.add_argument("--sample-read-size", help="Comma separated list of sample name and size (bp). If provided, sample read counting is skipped")
    coassemble_parser.add_argument("--genome-transcripts", nargs='+', help="Genome transcripts for reference database")
    coassemble_parser.add_argument("--genome-transcripts-list", help="Genome transcripts for reference database newline separated")
    coassemble_parser.add_argument("--genome-singlem", help="Combined SingleM otu tables for genome transcripts. If provided, genome SingleM is skipped")
    # Clustering options
    coassemble_parser.add_argument("--taxa-of-interest", help="Only consider sequences from this GTDB taxa (e.g. p__Planctomycetota) [default: all]")
    coassemble_parser.add_argument("--appraise-sequence-identity", type=int, help="Minimum sequence identity for SingleM appraise against reference database [default: 89%]", default=0.89)
    coassemble_parser.add_argument("--min-sequence-coverage", type=int, help="Minimum combined coverage for sequence inclusion [default: 10]", default=10)
    coassemble_parser.add_argument("--single-assembly", action="store_true", help="Skip appraise to discover samples to differential abundance binning. Forces --num-coassembly-samples and --max-coassembly-samples to 1")
    coassemble_parser.add_argument("--num-coassembly-samples", type=int, help="Number of samples per coassembly cluster [default: 2]", default=2)
    coassemble_parser.add_argument("--max-coassembly-samples", type=int, help="Upper bound for number of samples per coassembly cluster [default: --num-coassembly-samples]", default=None)
    coassemble_parser.add_argument("--max-coassembly-size", type=int, help="Maximum size (Gbp) of coassembly cluster [default: None]", default=None)
    coassemble_parser.add_argument("--max-recovery-samples", type=int, help="Upper bound for number of related samples to use for differential abundance binning [default: 20]", default=20)
    # Coassembly options
    coassemble_parser.add_argument("--assemble-unmapped", action="store_true", help="Only assemble reads that do not map to reference genomes")
    coassemble_parser.add_argument("--abstract-options", action="store_true", help="Print Aviary commands with bash variables for OUTPUT_DIR, CPUS and MEMORY [default: hardcode arguments]")
    coassemble_parser.add_argument("--aviary-cores", type=int, help="Maximum number of cores for Aviary to use", default=16)
    coassemble_parser.add_argument("--aviary-memory", type=int, help="Maximum amount of memory for Aviary to use (Gigabytes)", default=250)
    # General options
    coassemble_parser.add_argument("--output", help="Output directory [default: .]", default="./")
    coassemble_parser.add_argument("--conda-prefix", help="Path to conda environment install location", default=None)
    coassemble_parser.add_argument("--cores", type=int, help="Maximum number of cores to use", default=1)
    coassemble_parser.add_argument("--dryrun", action="store_true", help="dry run workflow")
    coassemble_parser.add_argument("--snakemake-args", help="Additional commands to be supplied to snakemake in the form of a space-prefixed single string e.g. \" --quiet\"", default="")

    ###########################################################################

    evaluate_parser = main_parser.new_subparser("evaluate", "Evaluate coassembled bins")
    evaluate_parser.add_argument("--cluster-output", help="Output dir from cluster subcommand", required=True)
    evaluate_parser.add_argument("--aviary-outputs", nargs='+', help="Output dir from Aviary coassembly and recover commands produced by coassemble subcommand", required=True)
    evaluate_parser.add_argument("--singlem-metapackage", help="SingleM metapackage for sequence searching")
    evaluate_parser.add_argument("--output", help="Output directory [default: .]", default="./")
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

    if args.subparser_name == "coassemble":
        coassemble(args)
    elif args.subparser_name == "evaluate":
        evaluate(args)

if __name__ == "__main__":
    sys.exit(main())
