##########################
### aviary_commands.py ###
##########################
# Author: Samuel Aroney

import pandas as pd

def produce_command(row, reads_1, reads_2, output_dir, threads, memory):
    coassembly_sample_names = row["samples"]
    coassembly_samples_1 = [reads_1[sample] for sample in coassembly_sample_names]
    coassembly_samples_2= [reads_2[sample] for sample in coassembly_sample_names]

    all_samples_names = row["recover_samples"]
    all_samples_1 = [reads_1[sample] for sample in all_samples_names]
    all_samples_2 = [reads_2[sample] for sample in all_samples_names]

    coassembly = row["coassembly"]
    output_dir = output_dir + "/coassemble"

    aviary_assemble = ("aviary assemble "
    f"-1 {' '.join(coassembly_samples_1)} "
    f"-2 {' '.join(coassembly_samples_2)} "
    f"--output {output_dir}/{coassembly}/assemble "
    f"-n {threads} "
    f"-t {threads} "
    f"-m {memory} "
    f"&> {output_dir}/logs/{coassembly}_assemble.log ")

    aviary_recover = ("aviary recover "
    f"--assembly {output_dir}/{coassembly}/assemble/assembly/final_contigs.fasta "
    f"-1 {' '.join(all_samples_1)} "
    f"-2 {' '.join(all_samples_2)} "
    f"--output {output_dir}/{coassembly}/recover "
    f"-n {threads} "
    f"-t {threads} "
    f"-m {memory} "
    f"&> {output_dir}/logs/{coassembly}_recover.log ")

    return pd.Series([aviary_assemble, aviary_recover])

def pipeline(coassemblies, reads_1, reads_2, output_dir, threads, memory):
    coassemblies["samples"] = coassemblies["samples"].apply(lambda item: item.split(","))
    coassemblies["recover_samples"] = coassemblies["recover_samples"].apply(lambda item: item.split(","))
    coassemblies["coassembly"] = coassemblies.reset_index()["index"].apply(lambda x: "coassembly_" + str(x))
    coassemblies[["assemble", "recover"]] = coassemblies.apply(
        produce_command,
        axis=1,
        reads_1=reads_1, reads_2=reads_2, output_dir=output_dir, threads=threads, memory=memory)

    return coassemblies

if __name__ == "__main__":
    coassemblies = pd.read_csv(snakemake.input.elusive_clusters, sep="\t")
    if len(coassemblies.index) == 0:
        with open(snakemake.output.coassemble_commands, "w") as f:
            pass
        with open(snakemake.output.recover_commands, "w") as f:
            pass
        print("No coassemblies to perform")
        exit(0)

    coassemblies = pipeline(
        coassemblies,
        snakemake.params.reads_1,
        snakemake.params.reads_2,
        snakemake.params.dir,
        snakemake.params.threads,
        snakemake.params.memory,
    )

    coassemblies["assemble"].to_csv(snakemake.output.coassemble_commands, sep="\t", index=False, header=False)
    coassemblies["recover"].to_csv(snakemake.output.recover_commands, sep="\t", index=False, header=False)
