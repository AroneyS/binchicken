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
        dir=directory(output_dir + "/mapping/{read}_coverm"),
        unmapped_bam=output_dir + "/mapping/{read}_unmapped.bam",
        reads_1=output_dir + "/mapping/{read}_unmapped.1.fq.gz",
        reads_2=output_dir + "/mapping/{read}_unmapped.2.fq.gz",
    log:
        logs_dir + "/mapping/{read}.log",
    params:
        genomes="{read}_reference.fna",
        reads_1=lambda wildcards: os.path.basename(config["reads_1"][wildcards.read]),
        coverm_log=logs_dir + "/mapping/{read}_coverm.log",
        view_log=logs_dir + "/mapping/{read}_view.log",
        fastq_log=logs_dir + "/mapping/{read}_fastq.log",
    threads:
        8
    conda:
        "env/coverm.yml"
    shell:
         "coverm make "
         "-r {input.genomes} "
         "-1 {input.reads_1} "
         "-2 {input.reads_2} "
         "-o {output.dir} "
         "-t {threads} "
         "&> {params.coverm_log} "
         "&& "
         "samtools view "
         "-@ $(({threads} - 1)) "
         "-b -f12 {output.dir}/{params.genomes}.{params.reads_1}.bam "
         "2> {params.view_log} "
         "> {output.unmapped_bam} "
         "&& "
         "samtools fastq "
         "-@ $(({threads} - 1)) "
         "{output.unmapped_bam} "
         "-1 {output.reads_1} "
         "-2 {output.reads_2} "
         "-0 /dev/null "
         "-s /dev/null "
         "-n "
         "&> {params.fastq_log} "

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
