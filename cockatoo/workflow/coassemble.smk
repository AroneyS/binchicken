#############
### Setup ###
#############
import os

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

rule all:
    input:
        output_dir + "/commands/coassemble_commands.sh",
        output_dir + "/commands/recover_commands.sh"

##################################
### Map reads to matching bins ###
##################################
rule collect_bins:
    input:
        appraise_binned=lambda wildcards: config["appraise_binned"][wildcards.read],
    output:
        output_dir + "/mapping/{read}_reference.fna",
    params:
        genomes=config["genomes"],
    script:
        "scripts/collect_reference_bins.py"

rule map_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read],
        genomes=output_dir + "/mapping/{read}_reference.fna",
    output:
        reads_1=output_dir + "/mapping/{read}_unmapped.1.fq.gz",
        reads_2=output_dir + "/mapping/{read}_unmapped.2.fq.gz",
    log:
        logs_dir + "/mapping/{read}.log"
    conda:
        "env/coverm.yml"
    shell:
        "touch {output.reads_1} {output.reads_2} "
        # "coverm make "
        # "-r {input.genomes} "
        # "--forward {input.reads_1} "
        # "--reverse {input.reads_2} "
        # "--genomes {input.genomes} "
        # "--mapped-forward {output.reads_1} "
        # "--mapped-reverse {output.reads_2} "
        # "&> {log}"

rule finish_mapping:
    input:
        mapped_reads_1.values(),
        mapped_reads_2.values(),
    output:
        output_dir + "/mapping/done"
    shell:
        "touch {output}"

##############################
### Create Aviary commands ###
##############################
rule coassemble_commands:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
        elusive_clusters=config["elusive_clusters"],
    output:
        coassemble_commands=output_dir + "/commands/coassemble_commands.sh",
        recover_commands=output_dir + "/commands/recover_commands.sh"
    params:
        reads_1 = mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"],
        reads_2 = mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"],
    log:
        logs_dir + "/coassemble_commands.log"
    script:
        "scripts/coassemble_commands.py"
