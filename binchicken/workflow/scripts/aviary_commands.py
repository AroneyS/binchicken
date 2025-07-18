#!/usr/bin/env python3
##########################
### aviary_commands.py ###
##########################
# Author: Samuel Aroney

import os
import polars as pl
from binchicken.binchicken import FAST_AVIARY_MODE
import argparse

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/aviary_commands.py \
#   --elusive-clusters {input.elusive_clusters} \
#   --coassemble-commands {output.coassemble_commands} \
#   --recover-commands {output.recover_commands} \
#   --reads-1 {input.reads_1} \
#   --reads-2 {input.reads_2} \
#   --dir {params.dir} \
#   --assemble-threads {params.assemble_threads} \
#   --assemble-memory {params.assemble_memory} \
#   --recover-threads {params.recover_threads} \
#   --recover-memory {params.recover_memory} \
#   --speed {params.speed} \
#   --threads {threads} \
#   --log {log}
# """

def pipeline(coassemblies, reads_1, reads_2, output_dir, assemble_threads, assemble_memory, recover_threads, recover_memory, fast=False):
    output = (
        coassemblies
        .with_columns(
            pl.col("samples").str.split(","),
            pl.col("recover_samples").str.split(","),
            )
        .with_columns(
            coassembly_samples_1 = pl.col("samples").list.eval(pl.element().replace(reads_1)),
            coassembly_samples_2 = pl.col("samples").list.eval(pl.element().replace(reads_2)),
            recover_samples_1 = pl.col("recover_samples").list.eval(pl.element().replace(reads_1)),
            recover_samples_2 = pl.col("recover_samples").list.eval(pl.element().replace(reads_2)),
            )
        .with_columns(
            assemble = pl.concat_str(
                pl.lit("aviary assemble --coassemble -1 "),
                pl.col("coassembly_samples_1").list.join(" "),
                pl.lit(" -2 "),
                pl.col("coassembly_samples_2").list.join(" "),
                pl.lit(" --output "),
                pl.lit(output_dir),
                pl.lit("/coassemble/"),
                pl.col("coassembly"),
                pl.lit("/assemble -n "),
                pl.lit(assemble_threads),
                pl.lit(" -t "),
                pl.lit(assemble_threads),
                pl.lit(" -m "),
                pl.lit(assemble_memory),
                pl.lit(" --skip-qc &> "),
                pl.lit(output_dir),
                pl.lit("/coassemble/logs/"),
                pl.col("coassembly"),
                pl.lit("_assemble.log ")
                ),
            recover = pl.concat_str(
                pl.lit("aviary recover --assembly "),
                pl.lit(output_dir),
                pl.lit("/coassemble/"),
                pl.col("coassembly"),
                pl.lit("/assemble/assembly/final_contigs.fasta -1 "),
                pl.col("recover_samples_1").list.join(" "),
                pl.lit(" -2 "),
                pl.col("recover_samples_2").list.join(" "),
                pl.lit(" --output "),
                pl.lit(output_dir),
                pl.lit("/coassemble/"),
                pl.col("coassembly"),
                pl.lit("/recover"),
                pl.when(pl.lit(fast)).then(pl.lit(" --binning-only --refinery-max-iterations 0")).otherwise(pl.lit("")),
                pl.lit(" -n "),
                pl.lit(recover_threads),
                pl.lit(" -t "),
                pl.lit(recover_threads),
                pl.lit(" -m "),
                pl.lit(recover_memory),
                pl.lit(" --skip-qc &> "),
                pl.lit(output_dir),
                pl.lit("/coassemble/logs/"),
                pl.col("coassembly"),
                pl.lit("_recover.log ")
                ),
            )
    )

    return output.select("assemble", "recover")

def main():
    parser = argparse.ArgumentParser(description="Aviary commands pipeline script.")
    parser.add_argument("--elusive-clusters", required=True, help="Path to elusive clusters input file")
    parser.add_argument("--coassemble-commands", required=True, help="Path to output coassemble commands file")
    parser.add_argument("--recover-commands", required=True, help="Path to output recover commands file")
    parser.add_argument("--reads-1", required=True, help="Named list file of read1")
    parser.add_argument("--reads-2", required=True, help="Named list file of read2")
    parser.add_argument("--dir", required=True, help="Output directory")
    parser.add_argument("--assemble-threads", type=int, required=True, help="Threads for assembly")
    parser.add_argument("--assemble-memory", required=True, help="Memory for assembly")
    parser.add_argument("--recover-threads", type=int, required=True, help="Threads for recovery")
    parser.add_argument("--recover-memory", required=True, help="Memory for recovery")
    parser.add_argument("--speed", default=None, help="Speed mode (FAST_AVIARY_MODE)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if args.log:
        import logging
        logging.basicConfig(
            filename=args.log,
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    else:
        import logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    coassemblies = pl.read_csv(args.elusive_clusters, separator="\t")
    if coassemblies.height == 0:
        with open(args.coassemble_commands, "w") as f:
            pass
        with open(args.recover_commands, "w") as f:
            pass
        print("No coassemblies to perform")
        exit(0)

    fast = args.speed == FAST_AVIARY_MODE

    reads_1 = {}
    with open(args.reads_1, "r") as f:
        for line in f:
            sample, read1 = line.strip().split("\t")
            reads_1[sample] = read1
    reads_2 = {}
    with open(args.reads_2, "r") as f:
        for line in f:
            sample, read2 = line.strip().split("\t")
            reads_2[sample] = read2

    coassemblies = pipeline(
        coassemblies,
        reads_1=reads_1,
        reads_2=reads_2,
        output_dir=args.dir,
        assemble_threads=args.assemble_threads,
        assemble_memory=args.assemble_memory,
        recover_threads=args.recover_threads,
        recover_memory=args.recover_memory,
        fast=fast,
    )

    coassemblies.select("assemble").write_csv(args.coassemble_commands, separator="\t", include_header=False)
    coassemblies.select("recover").write_csv(args.recover_commands, separator="\t", include_header=False)

if __name__ == "__main__":
    main()
