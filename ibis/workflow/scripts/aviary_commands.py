##########################
### aviary_commands.py ###
##########################
# Author: Samuel Aroney

import polars as pl
from ibis.ibis import FAST_AVIARY_MODE

def pipeline(coassemblies, reads_1, reads_2, output_dir, threads, memory, fast=False):
    output = (
        coassemblies
        .with_columns(
            pl.col("samples").str.split(","),
            pl.col("recover_samples").str.split(","),
            )
        .with_columns(
            coassembly_samples_1 = pl.col("samples").list.eval(pl.element().map_dict(reads_1)),
            coassembly_samples_2 = pl.col("samples").list.eval(pl.element().map_dict(reads_2)),
            recover_samples_1 = pl.col("recover_samples").list.eval(pl.element().map_dict(reads_1)),
            recover_samples_2 = pl.col("recover_samples").list.eval(pl.element().map_dict(reads_2)),
            )
        .with_columns(
            assemble = pl.concat_str(
                pl.lit("aviary assemble -1 "),
                pl.col("coassembly_samples_1").list.join(" "),
                pl.lit(" -2 "),
                pl.col("coassembly_samples_2").list.join(" "),
                pl.lit(" --output "),
                pl.lit(output_dir),
                pl.lit("/coassemble/"),
                pl.col("coassembly"),
                pl.lit("/assemble -n "),
                pl.lit(threads),
                pl.lit(" -t "),
                pl.lit(threads),
                pl.lit(" -m "),
                pl.lit(memory),
                pl.lit(" &> "),
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
                pl.when(pl.lit(fast)).then(pl.lit(" --workflow recover_mags_no_singlem --skip-binners maxbin concoct rosella --skip-abundances --refinery-max-iterations 0")).otherwise(pl.lit("")),
                pl.lit(" -n "),
                pl.lit(threads//2),
                pl.lit(" -t "),
                pl.lit(threads//2),
                pl.lit(" -m "),
                pl.lit(memory//2),
                pl.lit(" &> "),
                pl.lit(output_dir),
                pl.lit("/coassemble/logs/"),
                pl.col("coassembly"),
                pl.lit("_recover.log ")
                ),
            )
    )

    return output.select("assemble", "recover")

if __name__ == "__main__":
    coassemblies = pl.read_csv(snakemake.input.elusive_clusters, separator="\t")
    if coassemblies.height == 0:
        with open(snakemake.output.coassemble_commands, "w") as f:
            pass
        with open(snakemake.output.recover_commands, "w") as f:
            pass
        print("No coassemblies to perform")
        exit(0)

    fast = snakemake.params.speed == FAST_AVIARY_MODE

    coassemblies = pipeline(
        coassemblies,
        reads_1=snakemake.params.reads_1,
        reads_2=snakemake.params.reads_2,
        output_dir=snakemake.params.dir,
        threads=snakemake.params.threads,
        memory=snakemake.params.memory,
        fast=fast,
    )

    coassemblies.select("assemble").write_csv(snakemake.output.coassemble_commands, separator="\t", has_header=False)
    coassemblies.select("recover").write_csv(snakemake.output.recover_commands, separator="\t", has_header=False)
