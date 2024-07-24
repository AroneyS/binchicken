##########################
### aviary_commands.py ###
##########################
# Author: Samuel Aroney

import os
import polars as pl
from binchicken.binchicken import FAST_AVIARY_MODE

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

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

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
        assemble_threads=snakemake.params.assemble_threads,
        assemble_memory=snakemake.params.assemble_memory,
        recover_threads=snakemake.params.recover_threads,
        recover_memory=snakemake.params.recover_memory,
        fast=fast,
    )

    coassemblies.select("assemble").write_csv(snakemake.output.coassemble_commands, separator="\t", include_header=False)
    coassemblies.select("recover").write_csv(snakemake.output.recover_commands, separator="\t", include_header=False)
