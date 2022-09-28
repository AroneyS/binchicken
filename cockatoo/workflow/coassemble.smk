#############
### Setup ###
#############
import os

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

wildcard_constraints:
    coassembly = "coassembly_\d+"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

def get_coassemblies():
    with open(config["elusive_clusters"]) as f:
        coassemblies = [line.strip().split("\t")[6] for line in f if not line.endswith("coassembly\n")]
    return coassemblies

def get_reads(coassembly, read_type, reads_dict):
    if read_type == "assemble":
        column = 0
    elif read_type == "recover":
        column = 5
    else:
        raise ValueError("read_type must be 'assemble' or 'recover'")
    with open(config["elusive_clusters"]) as f:
        cluster = [line.strip() for line in f if line.endswith(coassembly + "\n")][0]

    read_files = [reads_dict[read] for read in cluster.split("\t")[column].split(",")]
    return read_files

rule all:
    input:
        output_dir + "/commands/coassemble_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/commands/recover_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/coassemblies.done" if config["run_aviary"] else [],

##################################
### Map reads to matching bins ###
##################################
rule collect_bins:
    input:
        appraise_binned=config["appraise_binned"],
    output:
        output_dir + "/mapping/{read}_reference.fna",
    params:
        genomes=config["genomes"],
        sample="{read}",
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
    threads:
        64
    params:
        reads_1 = mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"],
        reads_2 = mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"],
        dir=output_dir,
        memory=config["aviary_memory"],
        threads=config["aviary_threads"],
    log:
        logs_dir + "/coassemble_commands.log"
    script:
        "scripts/coassemble_commands.py"

###########################
### Run Aviary commands ###
###########################
rule aviary_assemble:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
        reads_1=lambda wildcards: get_reads(wildcards.coassembly, 'assemble',
                                            mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"]
                                            ),
        reads_2=lambda wildcards: get_reads(wildcards.coassembly, 'assemble',
                                            mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"]
                                            ),
    output:
        dir=directory(output_dir + "/assemble/{coassembly}"),
        assembly=output_dir + "/assemble/{coassembly}/assembly/final_contigs.fasta",
    log:
        logs_dir + "/assemble/{coassembly}.log"
    threads:
        64
    params:
        memory=config["aviary_memory"],
    conda:
        "env/aviary.yml"
    shell:
        "aviary assemble "
        "-1 {input.reads_1} "
        "-2 {input.reads_2} "
        "--output {output.dir} "
        "-n {threads} "
        "-m {params.memory} "
        "&> {log} "

rule aviary_recover:
    input:
        assembly=output_dir + "/assemble/{coassembly}/assembly/final_contigs.fasta",
        reads_1=lambda wildcards: get_reads(wildcards.coassembly, 'recover',
                                            mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"]
                                            ),
        reads_2=lambda wildcards: get_reads(wildcards.coassembly, 'recover',
                                            mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"]
                                            ),
    output:
        directory(output_dir + "/recover/{coassembly}"),
    log:
        logs_dir + "/recover/{coassembly}.log"
    threads:
        64
    params:
        memory=config["aviary_memory"],
    conda:
        "env/aviary.yml"
    shell:
        "aviary recover "
        "--assembly {input.assembly} "
        "-1 {input.reads_1} "
        "-2 {input.reads_2} "
        "--output {output} "
        "-n {threads} "
        "-m {params.memory} "
        "&> {log} "

rule collate_coassemblies:
    input:
        expand(output_dir + "/recover/{coassembly}", coassembly=get_coassemblies()),
    output:
        output_dir + "/coassemblies.done"
    shell:
        "touch {output}"
