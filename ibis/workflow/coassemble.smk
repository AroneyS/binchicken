#############
### Setup ###
#############
ruleorder: no_genomes > query_processing > singlem_appraise

import os
import pandas as pd

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

def get_reads(wildcards, forward=True, version=None):
    version = version if version else wildcards.version
    if version == "":
        if forward:
            return config["reads_1"]
        else:
            return config["reads_2"]
    elif version == "unmapped_":
        if forward:
            return mapped_reads_1
        else:
            return mapped_reads_2
    else:
        raise ValueError("Version should be empty or unmapped")

def get_cat(wildcards):
    return "zcat" if [r for r in get_reads(wildcards).values()][0].endswith(".gz") else "cat"

def get_reads_coassembly(wildcards, forward=True, recover=False):
    checkpoint_output = checkpoints.cluster_graph.get(**wildcards).output[0]
    elusive_clusters = pd.read_csv(checkpoint_output, sep = "\t")

    if recover:
        sample_names = elusive_clusters[elusive_clusters["coassembly"] == wildcards.coassembly]["recover_samples"].iloc[0]
    else:
        sample_names = elusive_clusters[elusive_clusters["coassembly"] == wildcards.coassembly]["samples"].iloc[0]

    sample_names = sample_names.split(",")

    if config["assemble_unmapped"]:
        version = "unmapped_"
    else:
        version = ""

    reads = get_reads(wildcards, forward=forward, version=version)
    return [reads[n] for n in sample_names]

def get_coassemblies(wildcards):
    checkpoint_output = checkpoints.cluster_graph.get().output[0]
    elusive_clusters = pd.read_csv(checkpoint_output, sep = "\t")
    coassemblies = elusive_clusters["coassembly"].to_list()

    return [output_dir + f"/coassemble/{c}/recover" for c in coassemblies]

rule all:
    input:
        output_dir + "/target/elusive_clusters.tsv",
        output_dir + "/commands/coassemble_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/commands/recover_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/commands/done" if config["run_aviary"] else [],
        output_dir + "/summary.tsv",

rule summary:
    input:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
        read_size=output_dir + "/unmapped_read_size.csv" if config["assemble_unmapped"] else [],
    output:
        summary=output_dir + "/summary.tsv",
    script:
        "scripts/summarise_coassemblies.py"

#####################
### SingleM reads ###
#####################
rule singlem_pipe_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read],
    output:
        output_dir + "/pipe/{read}_read.otu_table.tsv"
    log:
        logs_dir + "/pipe/{read}_read.log"
    params:
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem pipe "
        "--forward {input.reads_1} "
        "--reverse {input.reads_2} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

#######################
### SingleM genomes ###
#######################
rule genome_transcripts:
    input:
        lambda wildcards: config["genomes"][wildcards.genome],
    output:
        output_dir + "/transcripts/{genome}_protein.fna"
    log:
        logs_dir + "/transcripts/{genome}_protein.log"
    params:
        prodigal_meta = "-p meta" if config["prodigal_meta"] else ""
    conda:
        "env/prodigal.yml"
    shell:
        "prodigal "
        "-i {input} "
        "-d {output} "
        "{params.prodigal_meta} "
        "&> {log} "

rule singlem_pipe_genomes:
    input:
        output_dir + "/transcripts/{genome}_protein.fna"
    output:
        output_dir + "/pipe/{genome}_bin.otu_table.tsv"
    log:
        logs_dir + "/pipe/{genome}_bin.log"
    params:
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem pipe "
        "--forward {input} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

rule singlem_summarise_genomes:
    input:
        expand(output_dir + "/pipe/{genome}_bin.otu_table.tsv", genome=config["genomes"])
    output:
        output_dir + "/summarise/bins_summarised.otu_table.tsv"
    log:
        logs_dir + "/summarise/genomes.log"
    params:
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--exclude-off-target-hits "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

########################
### SingleM appraise ###
########################
rule singlem_appraise:
    input:
        reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        bins=output_dir + "/summarise/bins_summarised.otu_table.tsv",
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/binned.otu_table.tsv",
    log:
        logs_dir + "/appraise/appraise.log"
    params:
        sequence_identity=config["appraise_sequence_identity"],
        singlem_metapackage=config["singlem_metapackage"],
    conda:
        "env/singlem.yml"
    shell:
        "singlem appraise "
        "--metagenome-otu-tables {input.reads} "
        "--genome-otu-tables {input.bins} "
        "--metapackage {params.singlem_metapackage} "
        "--output-unaccounted-for-otu-table {output.unbinned} "
        "--output-binned-otu-table {output.binned} "
        "--imperfect "
        "--sequence-identity {params.sequence_identity} "
        "--output-found-in "
        "&> {log}"

###################################
### SingleM query (alternative) ###
###################################
rule query_processing:
    input:
        pipe_reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        query_reads=expand(output_dir + "/query/{read}_query.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/binned.otu_table.tsv",
    log:
        logs_dir + "/query/processing.log"
    params:
        sequence_identity=config["appraise_sequence_identity"],
        window_size=60,
    threads:
        64
    script:
        "scripts/query_processing.py"

################################
### No genomes (alternative) ###
################################
rule no_genomes:
    input:
        reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv" if config["no_genomes"] else [],
        binned=output_dir + "/appraise/binned.otu_table.tsv" if config["no_genomes"] else [],
    log:
        logs_dir + "/appraise/appraise.log"
    script:
        "scripts/no_genomes.py"

######################
### Target elusive ###
######################
rule count_bp_reads:
    input:
        reads_1 = lambda wildcards: get_reads(wildcards).values(),
        reads_2 = lambda wildcards: get_reads(wildcards, forward = False).values()
    output:
        output_dir + "/{version,.*}read_size.csv"
    params:
        names=list(config["reads_1"].keys()),
        cat=get_cat,
    threads:
        8
    shell:
        "parallel -k -j {threads} "
        "echo -n {{1}}, '&&' "
        "{params.cat} {{2}} {{3}} '|' sed -n 2~4p '|' tr -d '\"\n\"' '|' wc -m "
        "::: {params.names} :::+ {input.reads_1} :::+ {input.reads_2} "
        "> {output}"

rule target_elusive:
    input:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv"
    output:
        output_edges=output_dir + "/target/elusive_edges.tsv",
        output_targets=output_dir + "/target/targets.tsv",
    log:
        logs_dir + "/elusive_targets.log"
    params:
        min_coassembly_coverage=config["min_coassembly_coverage"],
        taxa_of_interest=config["taxa_of_interest"],
    threads:
        64
    script:
        "scripts/target_elusive.py"

checkpoint cluster_graph:
    input:
        elusive_edges=output_dir + "/target/elusive_edges.tsv",
        read_size=output_dir + "/read_size.csv",
    output:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv"
    log:
        logs_dir + "/cluster_graph.log"
    params:
        max_coassembly_size=config["max_coassembly_size"],
        num_coassembly_samples=config["num_coassembly_samples"],
        max_coassembly_samples=config["max_coassembly_samples"],
        max_recovery_samples=config["max_recovery_samples"],
    threads:
        64
    script:
        "scripts/cluster_graph.py"

#######################
### SRA downloading ###
#######################
rule download_sra:
    output:
        directory(output_dir + "/sra")
    threads:
        64
    params:
        sra=" ".join(config["sra"]) if config["sra"] else ""
    conda:
        "env/kingfisher.yml"
    log:
        logs_dir + "/sra/kingfisher.log"
    shell:
        "mkdir -p {output} && "
        "cd {output} && "
        "kingfisher get "
        "-r {params.sra} "
        "-f fastq.gz "
        "-m ena-ascp ena-ftp prefetch aws-http aws-cp "
        "-t {threads} "
        "&> {log} "

#####################################
### Map reads to matching genomes ###
#####################################
rule collect_genomes:
    input:
        appraise_binned=output_dir + "/appraise/binned.otu_table.tsv",
        appraise_unbinned=output_dir + "/appraise/unbinned.otu_table.tsv",
    output:
        temp(output_dir + "/mapping/{read}_reference.fna"),
    params:
        genomes=config["genomes"],
        sample="{read}",
        min_appraised=config["unmapping_min_appraised"],
    script:
        "scripts/collect_reference_bins.py"

rule map_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read],
        genomes=output_dir + "/mapping/{read}_reference.fna",
    output:
        dir=temp(directory(output_dir + "/mapping/{read}_coverm")),
    log:
        logs_dir + "/mapping/{read}_coverm.log",
    threads:
        32
    conda:
        "env/coverm.yml"
    shell:
        "coverm make "
        "-r {input.genomes} "
        "-1 {input.reads_1} "
        "-2 {input.reads_2} "
        "-o {output.dir} "
        "-t {threads} "
        "&> {log} "

rule filter_bam_files:
    input:
        output_dir + "/mapping/{read}_coverm",
    output:
        temp(output_dir + "/mapping/{read}_unmapped.bam"),
    log:
        logs_dir + "/mapping/{read}_filter.log",
    params:
        genomes="{read}_reference.fna",
        reads_1=lambda wildcards: os.path.basename(config["reads_1"][wildcards.read]),
        sequence_identity=config["unmapping_max_identity"],
    threads:
        32
    conda:
        "env/coverm.yml"
    shell:
        "coverm filter "
        "-b {input}/{params.genomes}.{params.reads_1}.bam "
        "-o {output} "
        "--inverse "
        "--min-read-percent-identity {params.sequence_identity} "
        "&> {log}"

rule bam_to_fastq:
    input:
        output_dir + "/mapping/{read}_unmapped.bam",
    output:
        reads_1=output_dir + "/mapping/{read}_unmapped.1.fq.gz",
        reads_2=output_dir + "/mapping/{read}_unmapped.2.fq.gz",
    log:
        logs_dir + "/mapping/{read}_fastq.log",
    threads:
        32
    conda:
        "env/coverm.yml"
    shell:
        "samtools fastq "
        "-@ $(({threads} - 1)) "
        "{input} "
        "-1 {output.reads_1} "
        "-2 {output.reads_2} "
        "-0 /dev/null "
        "-s /dev/null "
        "-n "
        "&> {log} "

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
rule aviary_commands:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
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
        logs_dir + "/aviary_commands.log"
    script:
        "scripts/aviary_commands.py"

#########################################
### Run Aviary commands (alternative) ###
#########################################
rule aviary_assemble:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
    output:
        dir=directory(output_dir + "/coassemble/{coassembly}/assemble"),
        assembly=output_dir + "/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta",
    params:
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False),
        memory=config["aviary_memory"],
    threads:
        threads=config["aviary_threads"]
    log:
        logs_dir + "/aviary/{coassembly}_assemble.log"
    conda:
        config["aviary_conda"]
    shell:
        "aviary assemble "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {output.dir} "
        "-n {threads} "
        "-t {threads} "
        "-m {params.memory} "
        "&> {log} "

rule aviary_recover:
    input:
        assembly=output_dir + "/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta",
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
    output:
        directory(output_dir + "/coassemble/{coassembly}/recover")
    params:
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards, recover=True),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False, recover=True),
        memory=config["aviary_memory"],
        skip_binners = "--skip-binners rosella semibin metabat1 vamb concoct maxbin2 maxbin" if config["test"] else "",
    threads:
        config["aviary_threads"]
    log:
        logs_dir + "/aviary/{coassembly}_recover.log"
    conda:
        config["aviary_conda"]
    shell:
        "aviary recover "
        "--assembly {input.assembly} "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {output} "
        "{params.skip_binners} "
        "-n {threads} "
        "-t {threads} "
        "-m {params.memory} "
        "&> {log} "

rule aviary_combine:
    input:
        get_coassemblies,
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
        mapping=output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
    output:
        output_dir + "/commands/done",
    shell:
        "touch {output} "
