#############
### Setup ###
#############
ruleorder: no_genomes > query_processing > update_appraise > singlem_appraise
ruleorder: mock_download_sra > download_sra
localrules: all, summary, singlem_summarise_genomes, singlem_appraise_filtered, no_genomes, mock_download_sra, compile_sra_qc, collect_genomes, finish_mapping, aviary_commands, aviary_combine

import os
import pandas as pd
from ibis.ibis import FAST_AVIARY_MODE

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

def get_genomes(wildcards, version=None):
    version = version if version else wildcards.version
    if version == "":
        return expand(output_dir + "/pipe/{genome}_bin.otu_table.tsv", genome=config["genomes"])
    elif version == "new_":
        return expand(output_dir + "/pipe/{genome}_bin.otu_table.tsv", genome=config["new_genomes"])
    else:
        raise ValueError("Version should be empty or 'new'")

def get_reads(wildcards, forward=True, version=None):
    version = version if version else wildcards.version
    if version == "" or version == "whole":
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
        raise ValueError("Version should be empty, 'whole' or 'unmapped'")

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
        version = "whole"

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
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
        read_size = output_dir + "/unmapped_read_size.csv" if config["assemble_unmapped"] else [],
    output:
        summary = output_dir + "/summary.tsv",
    script:
        "scripts/summarise_coassemblies.py"

#####################
### SingleM reads ###
#####################
rule singlem_pipe_reads:
    input:
        reads_1 = lambda wildcards: config["reads_1"][wildcards.read],
        reads_2 = lambda wildcards: config["reads_2"][wildcards.read],
    output:
        output_dir + "/pipe/{read}_read.otu_table.tsv"
    log:
        logs_dir + "/pipe/{read}_read.log"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
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
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
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
        singlem_metapackage = config["singlem_metapackage"]
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
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
        lambda wildcards: get_genomes(wildcards)
    output:
        output_dir + "/summarise/{version,.*}bins_summarised.otu_table.tsv"
    log:
        logs_dir + "/summarise/{version,.*}genomes.log"
    params:
        singlem_metapackage = config["singlem_metapackage"]
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
        reads = expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        bins = output_dir + "/summarise/bins_summarised.otu_table.tsv",
    output:
        unbinned = temp(output_dir + "/appraise/unbinned_raw.otu_table.tsv"),
        binned = temp(output_dir + "/appraise/binned_raw.otu_table.tsv"),
    log:
        logs_dir + "/appraise/appraise.log"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        singlem_metapackage = config["singlem_metapackage"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
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

rule singlem_appraise_filtered:
    input:
        unbinned = output_dir + "/appraise/unbinned_raw.otu_table.tsv",
        binned = output_dir + "/appraise/binned_raw.otu_table.tsv",
    output:
        unbinned = output_dir + "/appraise/unbinned.otu_table.tsv",
        binned = output_dir + "/appraise/binned.otu_table.tsv",
    params:
        bad_package = "S3.18.EIF_2_alpha",
    shell:
        "grep -v {params.bad_package} {input.unbinned} > {output.unbinned} && "
        "grep -v {params.bad_package} {input.binned} > {output.binned}"

#####################################
### Update appraise (alternative) ###
#####################################
rule update_appraise:
    input:
        unbinned = output_dir + "/appraise/unbinned_prior.otu_table.tsv",
        binned = output_dir + "/appraise/binned_prior.otu_table.tsv",
        bins = output_dir + "/summarise/new_bins_summarised.otu_table.tsv",
    output:
        unbinned = temp(output_dir + "/appraise/unbinned_raw.otu_table.tsv"),
        binned = temp(output_dir + "/appraise/binned_raw.otu_table.tsv"),
    log:
        logs_dir + "/appraise/appraise.log"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        singlem_metapackage = config["singlem_metapackage"],
        new_binned = output_dir + "/appraise/binned_new.otu_table.tsv",
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
    conda:
        "env/singlem.yml"
    shell:
        "singlem appraise "
        "--metagenome-otu-tables {input.unbinned} "
        "--genome-otu-tables {input.bins} "
        "--metapackage {params.singlem_metapackage} "
        "--output-unaccounted-for-otu-table {output.unbinned} "
        "--output-binned-otu-table {params.new_binned} "
        "--imperfect "
        "--sequence-identity {params.sequence_identity} "
        "--output-found-in "
        "&> {log} "
        "&& cat {input.binned} <(tail -n+2 {params.new_binned}) > {output.binned} "

###################################
### SingleM query (alternative) ###
###################################
rule query_processing:
    input:
        pipe_reads = expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        query_reads = expand(output_dir + "/query/{read}_query.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned = temp(output_dir + "/appraise/unbinned_raw.otu_table.tsv"),
        binned = temp(output_dir + "/appraise/binned_raw.otu_table.tsv"),
    log:
        logs_dir + "/query/processing.log"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        window_size = 60,
        taxa_of_interest = config["taxa_of_interest"],
    threads: 1
    resources:
        mem_mb=8000,
        runtime = "24h",
    script:
        "scripts/query_processing.py"

################################
### No genomes (alternative) ###
################################
rule no_genomes:
    input:
        reads = expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned = temp(output_dir + "/appraise/unbinned_raw.otu_table.tsv") if config["no_genomes"] else [],
        binned = temp(output_dir + "/appraise/binned_raw.otu_table.tsv") if config["no_genomes"] else [],
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
        reads_2 = lambda wildcards: get_reads(wildcards, forward=False).values()
    output:
        output_dir + "/{version,.*}read_size.csv"
    params:
        names = list(config["reads_1"].keys()),
        cat = get_cat,
    threads: 8
    resources:
        mem_mb=64000,
        runtime = "24h",
    shell:
        "parallel -k -j {threads} "
        "echo -n {{1}}, '&&' "
        "{params.cat} {{2}} {{3}} '|' sed -n 2~4p '|' tr -d '\"\n\"' '|' wc -m "
        "::: {params.names} :::+ {input.reads_1} :::+ {input.reads_2} "
        "> {output}"

rule target_elusive:
    input:
        unbinned = output_dir + "/appraise/unbinned.otu_table.tsv"
    output:
        output_edges = output_dir + "/target/elusive_edges.tsv",
        output_targets = output_dir + "/target/targets.tsv",
    params:
        min_coassembly_coverage = config["min_coassembly_coverage"],
        max_coassembly_samples = config["max_coassembly_samples"],
        taxa_of_interest = config["taxa_of_interest"],
    threads: 32
    resources:
        mem_mb=250*1000,
        runtime = "24h",
    script:
        "scripts/target_elusive.py"

checkpoint cluster_graph:
    input:
        elusive_edges = output_dir + "/target/elusive_edges.tsv",
        read_size = output_dir + "/read_size.csv",
    output:
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv"
    params:
        max_coassembly_size = config["max_coassembly_size"],
        num_coassembly_samples = config["num_coassembly_samples"],
        max_coassembly_samples = config["max_coassembly_samples"],
        max_recovery_samples = config["max_recovery_samples"],
        exclude_coassemblies = config["exclude_coassemblies"],
    threads: 64
    resources:
        mem_mb=500*1000,
        runtime = "168h",
    script:
        "scripts/cluster_graph.py"

#######################
### SRA downloading ###
#######################
rule download_sra:
    output:
        directory(output_dir + "/sra")
    params:
        sra = " ".join(config["sra"]) if config["sra"] else ""
    threads: 32
    resources:
        mem_mb=250*1000,
        runtime = "48h",
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
        "-m ena-ftp prefetch ena-ascp aws-http aws-cp "
        "-t {threads} "
        "&> {log} "

rule mock_download_sra:
    output:
        directory(output_dir + "/sra") if config["mock_sra"] else []
    threads: 8
    params:
        sra_u = workflow.basedir + "/../../test/data/sra/" + config["sra"][0] + ".fastq.gz" if config["sra"] else "",
        sra_f = " ".join([workflow.basedir + "/../../test/data/sra/" + s + "_1.fastq.gz" for s in config["sra"][1:]]) if config["sra"] else "",
        sra_r = " ".join([workflow.basedir + "/../../test/data/sra/" + s + "_2.fastq.gz" for s in config["sra"][1:]]) if config["sra"] else "",
    conda:
        "env/kingfisher.yml"
    log:
        logs_dir + "/sra/kingfisher.log"
    shell:
        "mkdir -p {output} && "
        "cp {params.sra_u} {params.sra_f} {params.sra_r} {output}"

rule sra_qc:
    output:
        out1 = output_dir + "/sra_qc/{sra}_1.fastq.gz",
        out2 = output_dir + "/sra_qc/{sra}_2.fastq.gz",
    params:
        in1 = output_dir + "/sra/{sra}_1.fastq.gz",
        in2 = output_dir + "/sra/{sra}_2.fastq.gz",
        stats = output_dir + "/sra_qc/{sra}_qc.stats",
    threads: 16
    resources:
        mem_mb=125*1000,
        runtime = "24h",
    conda:
        "env/bbtools.yml"
    log:
        logs_dir + "/sra_qc/{sra}_qc.log"
    shell:
        "bbduk.sh "
        "in1={params.in1} in2={params.in2} "
        "out1={output.out1} out2={output.out2} "
        "ref=adapters,phix "
        "ktrim=r k=23 mink=11 hdist=1 "
        "qtrim=r trimq=10 "
        "minlen=30 "
        "stats={params.stats} statscolumns=3 "
        "t={threads} "
        "&> {log} "

rule compile_sra_qc:
    input:
        expand(output_dir + "/sra_qc/{sra}_1.fastq.gz", sra=config["sra"]),
        expand(output_dir + "/sra_qc/{sra}_2.fastq.gz", sra=config["sra"]),
    output:
        done = output_dir + "/sra_qc/done"
    shell:
        "touch {output.done} "

#####################################
### Map reads to matching genomes ###
#####################################
rule collect_genomes:
    input:
        appraise_binned = output_dir + "/appraise/binned.otu_table.tsv",
        appraise_unbinned = output_dir + "/appraise/unbinned.otu_table.tsv",
    output:
        temp(output_dir + "/mapping/{read}_reference.fna"),
    params:
        genomes = config["genomes"],
        sample = "{read}",
        min_appraised = config["unmapping_min_appraised"],
    script:
        "scripts/collect_reference_bins.py"

rule map_reads:
    input:
        reads_1 = lambda wildcards: config["reads_1"][wildcards.read],
        reads_2 = lambda wildcards: config["reads_2"][wildcards.read],
        genomes = output_dir + "/mapping/{read}_reference.fna",
    output:
        dir = temp(directory(output_dir + "/mapping/{read}_coverm")),
    threads: 16
    resources:
        mem_mb=125*1000,
        runtime = "24h",
    log:
        logs_dir + "/mapping/{read}_coverm.log",
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
    params:
        genomes = "{read}_reference.fna",
        reads_1 = lambda wildcards: os.path.basename(config["reads_1"][wildcards.read]),
        sequence_identity = config["unmapping_max_identity"],
        alignment_percent = config["unmapping_max_alignment"],
    threads: 16
    resources:
        mem_mb=125*1000,
        runtime = "24h",
    log:
        logs_dir + "/mapping/{read}_filter.log",
    conda:
        "env/coverm.yml"
    shell:
        "coverm filter "
        "-b {input}/{params.genomes}.{params.reads_1}.bam "
        "-o {output} "
        "--inverse "
        "--min-read-percent-identity-pair {params.sequence_identity} "
        "--min-read-aligned-percent-pair {params.alignment_percent} "
        "--proper-pairs-only "
        "&> {log}"

rule bam_to_fastq:
    input:
        output_dir + "/mapping/{read}_unmapped.bam",
    output:
        reads_1 = output_dir + "/mapping/{read}_unmapped.1.fq.gz",
        reads_2 = output_dir + "/mapping/{read}_unmapped.2.fq.gz",
    threads: 16
    resources:
        mem_mb=125*1000,
        runtime = "24h",
    log:
        logs_dir + "/mapping/{read}_fastq.log",
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
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
    output:
        coassemble_commands = output_dir + "/commands/coassemble_commands.sh",
        recover_commands = output_dir + "/commands/recover_commands.sh"
    threads: 8
    params:
        reads_1 = mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"],
        reads_2 = mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"],
        dir = output_dir,
        memory = config["aviary_memory"],
        threads = config["aviary_threads"],
        speed = config["aviary_speed"],
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
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
    output:
        dir = directory(output_dir + "/coassemble/{coassembly}/assemble"),
        assembly = output_dir + "/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta",
    params:
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False),
        dryrun = "--dryrun" if config["aviary_dryrun"] else "",
        drymkdir = "&& mkdir -p "+output_dir+"/coassemble/{coassembly}/assemble/assembly" if config["aviary_dryrun"] else "",
        drytouch = "&& touch "+output_dir+"/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta" if config["aviary_dryrun"] else "",
        conda_prefix = config["conda_prefix"] if config["conda_prefix"] else ".",
    threads:
        threads = config["aviary_threads"]
    resources:
        mem_mb = int(config["aviary_memory"]*1000),
        runtime = "96h",
        assembler = lambda wildcards, attempt: "" if attempt == 1 else "--use-megahit",
    log:
        logs_dir + "/aviary/{coassembly}_assemble.log"
    conda:
        "env/aviary.yml"
    shell:
        "GTDBTK_DATA_PATH=. "
        "CHECKM2DB=. "
        "EGGNOG_DATA_DIR=. "
        "CONDA_ENV_PATH={params.conda_prefix} "
        "aviary assemble "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {output.dir} "
        "-n {threads} "
        "-t {threads} "
        "-m {resources.mem_mb} "
        "{resources.assembler} "
        "{params.dryrun} "
        "&> {log} "
        "{params.drymkdir} "
        "{params.drytouch} "

rule aviary_recover:
    input:
        assembly = output_dir + "/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta",
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
    output:
        directory(output_dir + "/coassemble/{coassembly}/recover")
    params:
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards, recover=True),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False, recover=True),
        dryrun = "--dryrun" if config["aviary_dryrun"] else "",
        gtdbtk = config["aviary_gtdbtk"],
        checkm2 = config["aviary_checkm2"],
        conda_prefix = config["conda_prefix"] if config["conda_prefix"] else ".",
        fast = "--workflow recover_mags_no_singlem --skip-binners maxbin concoct rosella --skip-abundances --refinery-max-iterations 0" if config["aviary_speed"] == FAST_AVIARY_MODE else "",
    threads:
        int(config["aviary_threads"]/2)
    resources:
        mem_mb = int(config["aviary_memory"]*1000/2),
        runtime = "168h",
    log:
        logs_dir + "/aviary/{coassembly}_recover.log"
    conda:
        "env/aviary.yml"
    shell:
        "GTDBTK_DATA_PATH={params.gtdbtk} "
        "CHECKM2DB={params.checkm2} "
        "EGGNOG_DATA_DIR=. "
        "CONDA_ENV_PATH={params.conda_prefix} "
        "aviary recover "
        "--assembly {input.assembly} "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {output} "
        "{params.fast} "
        "-n {threads} "
        "-t {threads} "
        "-m {resources.mem_mb} "
        "{params.dryrun} "
        "&> {log} "

rule aviary_combine:
    input:
        get_coassemblies,
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
        mapping = output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
    output:
        output_dir + "/commands/done",
    shell:
        "touch {output} "
