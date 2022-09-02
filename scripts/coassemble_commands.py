#####################
### coassemble.py ###
#####################
# Author: Samuel Aroney
# Environment: coassembly
# Produce aviary assemble and recover commands for co-assembly.

import pandas as pd

#################
### Functions ###
#################
def produce_command(row):
    coassembly_sample_names = row["samples"]
    coassembly_samples_1 = [snakemake.config["reads_1"][sample] for sample in coassembly_sample_names]
    coassembly_samples_2= [snakemake.config["reads_2"][sample] for sample in coassembly_sample_names]

    all_samples_names = snakemake.config["reads_1"].keys()
    all_samples_1 = [snakemake.config["reads_1"][sample] for sample in all_samples_names]
    all_samples_2 = [snakemake.config["reads_2"][sample] for sample in all_samples_names]

    coassembly = row["coassembly"]

    aviary_assemble = ("aviary assemble "
    f"-1 {' '.join(coassembly_samples_1)} "
    f"-2 {' '.join(coassembly_samples_2)} "
    f"--output $OUTPUT_DIR/{coassembly}/assemble "
    f"-n $CPUS "
    f"-m $MEMORY "
    f"&> $OUTPUT_DIR/logs/{coassembly}_assemble.log")

    aviary_recover = ("aviary recover "
    f"--assembly $OUTPUT_DIR/{coassembly}/assemble/assembly/final_contigs.fasta "
    f"-1 {' '.join(all_samples_1)} "
    f"-2 {' '.join(all_samples_2)} "
    f"--output $OUTPUT_DIR/{coassembly}/recover "
    f"-n $CPUS "
    f"-m $MEMORY "
    f"&> $OUTPUT_DIR/logs/{coassembly}_recover.log ")

    return pd.Series([aviary_assemble, aviary_recover])

################
### Pipeline ###
################
coassemblies = pd.read_csv(snakemake.input.elusive_clusters, sep="\t")

# Produce coassembly_commands
coassemblies["samples"] = coassemblies["samples"].apply(lambda item: item.split(","))
coassemblies["coassembly"] = coassemblies.reset_index()["index"].apply(lambda x: "coassembly_" + str(x))
coassemblies[["assemble", "recover"]] = coassemblies.apply(produce_command, axis=1)

coassemblies["assemble"].to_csv(snakemake.output.coassemble_commands, sep="\t", index=False, header=False)
coassemblies["recover"].to_csv(snakemake.output.recover_commands, sep="\t", index=False, header=False)

