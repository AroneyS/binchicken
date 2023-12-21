#############
### Setup ###
#############
ruleorder: no_genomes > query_processing > update_appraise > singlem_appraise
ruleorder: mock_download_sra > download_sra

import os
import polars as pl
from binchicken.binchicken import FAST_AVIARY_MODE
os.umask(0o002)

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"
benchmarks_dir = output_dir + "/benchmarks"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

qc_reads_1 = {read: output_dir + f"/qc/{read}_1.fastq.gz" for read in config["reads_1"]}
qc_reads_2 = {read: output_dir + f"/qc/{read}_2.fastq.gz" for read in config["reads_2"]}

def get_mem_mb(wildcards, threads):
    return 8 * 1000 * threads

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
    elif version == "qc_":
        if forward:
            return qc_reads_1
        else:
            return qc_reads_2
    else:
        raise ValueError("Version should be empty, 'whole' or 'unmapped'")

def get_cat(wildcards):
    return "zcat" if [r for r in get_reads(wildcards).values()][0].endswith(".gz") else "cat"

def get_reads_coassembly(wildcards, forward=True, recover=False):
    checkpoint_output = checkpoints.cluster_graph.get(**wildcards).output[0]
    elusive_clusters = pl.read_csv(checkpoint_output, separator="\t")

    if recover:
        sample_names = elusive_clusters.filter(pl.col("coassembly") == wildcards.coassembly).get_column("recover_samples").to_list()[0]
    else:
        sample_names = elusive_clusters.filter(pl.col("coassembly") == wildcards.coassembly).get_column("samples").to_list()[0]

    sample_names = sample_names.split(",")

    if config["assemble_unmapped"]:
        version = "unmapped_"
    elif config["run_qc"]:
        version = "qc_"
    else:
        version = "whole"

    reads = get_reads(wildcards, forward=forward, version=version)
    return [reads[n] for n in sample_names]

def get_coassemblies(wildcards):
    checkpoint_output = checkpoints.cluster_graph.get().output[0]
    elusive_clusters = pl.read_csv(checkpoint_output, separator="\t")
    coassemblies = elusive_clusters.get_column("coassembly").to_list()

    return [output_dir + f"/coassemble/{c}/recover.done" for c in coassemblies]

rule all:
    input:
        output_dir + "/target/elusive_clusters.tsv",
        output_dir + "/commands/coassemble_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/commands/recover_commands.sh" if not config["run_aviary"] else [],
        output_dir + "/commands/done" if config["run_aviary"] else [],
        output_dir + "/summary.tsv",
    localrule: True

rule summary:
    input:
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
        read_size = output_dir + "/unmapped_read_size.csv" if config["assemble_unmapped"] else [],
    output:
        summary = output_dir + "/summary.tsv",
    localrule: True
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
    benchmark:
        benchmarks_dir + "/pipe/{read}_read.tsv"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    threads: 1
    resources:
        mem_mb=get_mem_mb,
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
    benchmark:
        benchmarks_dir + "/transcripts/{genome}_protein.tsv"
    params:
        prodigal_meta = "-p meta" if config["prodigal_meta"] else ""
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        runtime = "1h",
    group: "singlem_bins"
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
    benchmark:
        benchmarks_dir + "/pipe/{genome}_bin.tsv"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        runtime = "1h",
    group: "singlem_bins"
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
    benchmark:
        benchmarks_dir + "/summarise/{version,.*}genomes.tsv"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    localrule: True
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
    benchmark:
        benchmarks_dir + "/appraise/appraise.tsv"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        singlem_metapackage = config["singlem_metapackage"],
    threads: 1
    resources:
        mem_mb=get_mem_mb,
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
    localrule: True
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
    benchmark:
        benchmarks_dir + "/appraise/appraise.tsv"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        singlem_metapackage = config["singlem_metapackage"],
        new_binned = output_dir + "/appraise/binned_new.otu_table.tsv",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
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
    benchmark:
        benchmarks_dir + "/query/processing.tsv"
    params:
        sequence_identity = config["appraise_sequence_identity"],
        window_size = 60,
        taxa_of_interest = config["taxa_of_interest"],
    threads: 1
    resources:
        mem_mb=get_mem_mb,
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
    localrule: True
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
        mem_mb=get_mem_mb,
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
        samples = config["reads_1"],
    threads: 32
    resources:
        mem_mb=get_mem_mb,
        runtime = "24h",
    log:
        logs_dir + "/target/target_elusive.log"
    benchmark:
        benchmarks_dir + "/target/target_elusive.tsv"
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
        coassembly_samples = config["coassembly_samples"],
        exclude_coassemblies = config["exclude_coassemblies"],
    threads: 64
    resources:
        mem_mb=get_mem_mb,
        runtime = "168h",
    log:
        logs_dir + "/target/cluster_graph.log"
    benchmark:
        benchmarks_dir + "/target/cluster_graph.tsv"
    script:
        "scripts/cluster_graph.py"

#######################
### SRA downloading ###
#######################
rule download_read:
    output:
        done = output_dir + "/sra/{read}.done",
    params:
        dir = output_dir + "/sra",
        name = "{read}",
    threads: 4
    resources:
        mem_mb=get_mem_mb,
        runtime = "4h",
        downloading = 1,
    conda:
        "env/kingfisher.yml"
    log:
        logs_dir + "/sra/kingfisher_{read}.log"
    shell:
        "cd {params.dir} && "
        "rm -f {params.name}*.fastq.gz && "
        "kingfisher get "
        "-r {params.name} "
        "-f fastq.gz "
        "-m ena-ftp prefetch ena-ascp aws-http aws-cp "
        "-t {threads} "
        "&> {log} "
        "&& touch {output.done}"

rule download_sra:
    input:
        expand(output_dir + "/sra/{read}.done", read=config["sra"]) if config["sra"] else [],
    output:
        output_dir + "/sra/all_done"
    localrule: True
    shell:
        "touch {output}"

rule mock_download_sra:
    output:
        directory(output_dir + "/sra") if config["mock_sra"] else []
    threads: 8
    params:
        sra_u = workflow.basedir + "/../../test/data/sra/" + config["sra"][0] + ".fastq.gz" if config["sra"] else "",
        sra_f = " ".join([workflow.basedir + "/../../test/data/sra/" + s + "_1.fastq.gz" for s in config["sra"][1:]]) if config["sra"] else "",
        sra_r = " ".join([workflow.basedir + "/../../test/data/sra/" + s + "_2.fastq.gz" for s in config["sra"][1:]]) if config["sra"] else "",
    localrule: True
    conda:
        "env/kingfisher.yml"
    log:
        logs_dir + "/sra/kingfisher.log"
    shell:
        "mkdir -p {output} && "
        "cp {params.sra_u} {params.sra_f} {params.sra_r} {output}"

#####################################
### Map reads to matching genomes ###
#####################################
rule qc_reads:
    input:
        reads_1 = lambda wildcards: config["reads_1"][wildcards.read],
        reads_2 = lambda wildcards: config["reads_2"][wildcards.read],
    output:
        reads_1 = output_dir + "/qc/{read}_1.fastq.gz",
        reads_2 = output_dir + "/qc/{read}_2.fastq.gz",
        json = output_dir + "/qc/{read}.json",
        html = output_dir + "/qc/{read}.html",
    group: "unmapping"
    params:
        quality_cutoff = 15,
        unqualified_percent_limit = 40,
        min_length = 80,
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = "4h",
    log:
        logs_dir + "/mapping/{read}_qc.log"
    benchmark:
        benchmarks_dir + "/mapping/{read}_qc.tsv"
    conda:
        "env/fastp.yml"
    shell:
        "fastp "
        "-i {input.reads_1} "
        "-I {input.reads_2} "
        "-o {output.reads_1} "
        "-O {output.reads_2} "
        "-j {output.json} "
        "-h {output.html} "
        "-w {threads} "
        "-q {params.quality_cutoff} "
        "-u {params.unqualified_percent_limit} "
        "-l {params.min_length} "
        "&> {log}"

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
    localrule: True
    script:
        "scripts/collect_reference_bins.py"

rule map_reads:
    input:
        reads_1 = lambda wildcards: config["reads_1"][wildcards.read] if not config["run_qc"] else output_dir + "/qc/{read}_1.fastq.gz",
        reads_2 = lambda wildcards: config["reads_2"][wildcards.read] if not config["run_qc"] else output_dir + "/qc/{read}_2.fastq.gz",
        genomes = output_dir + "/mapping/{read}_reference.fna",
    output:
        dir = temp(directory(output_dir + "/mapping/{read}_coverm")),
    group: "unmapping"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = "12h",
    log:
        logs_dir + "/mapping/{read}_coverm.log",
    benchmark:
        benchmarks_dir + "/mapping/{read}_coverm.tsv"
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
    group: "unmapping"
    params:
        genomes = "{read}_reference.fna",
        reads_1 = lambda wildcards: os.path.basename(config["reads_1"][wildcards.read]),
        sequence_identity = config["unmapping_max_identity"],
        alignment_percent = config["unmapping_max_alignment"],
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = "4h",
    log:
        logs_dir + "/mapping/{read}_filter.log",
    benchmark:
        benchmarks_dir + "/mapping/{read}_filter.tsv"
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
    group: "unmapping"
    threads: 16
    resources:
        mem_mb=get_mem_mb,
        runtime = "4h",
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
    localrule: True
    shell:
        "touch {output}"

rule finish_qc:
    input:
        qc_reads_1.values(),
        qc_reads_2.values(),
    output:
        output_dir + "/qc/done"
    localrule: True
    shell:
        "touch {output}"

##############################
### Create Aviary commands ###
##############################
rule aviary_commands:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else output_dir + "/qc/done" if config["run_qc"] else [],
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
    output:
        coassemble_commands = output_dir + "/commands/coassemble_commands.sh",
        recover_commands = output_dir + "/commands/recover_commands.sh"
    threads: 8
    params:
        reads_1 = mapped_reads_1 if config["assemble_unmapped"] else qc_reads_1 if config["run_qc"] else config["reads_1"],
        reads_2 = mapped_reads_2 if config["assemble_unmapped"] else qc_reads_2 if config["run_qc"] else config["reads_2"],
        dir = output_dir,
        memory = config["aviary_memory"],
        threads = config["aviary_threads"],
        speed = config["aviary_speed"],
    localrule: True
    log:
        logs_dir + "/aviary_commands.log"
    script:
        "scripts/aviary_commands.py"

#########################################
### Run Aviary commands (alternative) ###
#########################################
rule aviary_assemble:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else output_dir + "/qc/done" if config["run_qc"] else [],
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
    output:
        dir = directory(output_dir + "/coassemble/{coassembly}/assemble"),
        assembly = output_dir + "/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta",
    params:
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False),
        dryrun = "--build" if config["build"] else "--dryrun" if config["aviary_dryrun"] else "",
        drymkdir = "&& mkdir -p "+output_dir+"/coassemble/{coassembly}/assemble/assembly" if config["aviary_dryrun"] else "",
        drytouch = "&& touch "+output_dir+"/coassemble/{coassembly}/assemble/assembly/final_contigs.fasta" if config["aviary_dryrun"] else "",
        conda_prefix = config["conda_prefix"] if config["conda_prefix"] else ".",
        tmpdir = config["tmpdir"],
    threads:
        threads = config["aviary_threads"]
    resources:
        mem_mb = int(config["aviary_memory"])*1000,
        mem_gb = int(config["aviary_memory"]),
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
        "SINGLEM_METAPACKAGE_PATH=. "
        "TMPDIR={params.tmpdir} "
        "aviary assemble "
        "--coassemble "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {output.dir} "
        "-n {threads} "
        "-t {threads} "
        "-m {resources.mem_gb} "
        "--skip-qc "
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
        output_dir + "/coassemble/{coassembly}/recover.done",
    params:
        output = output_dir + "/coassemble/{coassembly}/recover",
        reads_1 = lambda wildcards: get_reads_coassembly(wildcards, recover=True),
        reads_2 = lambda wildcards: get_reads_coassembly(wildcards, forward=False, recover=True),
        dryrun = "--build" if config["build"] else "--dryrun" if config["aviary_dryrun"] else "",
        gtdbtk = config["aviary_gtdbtk"],
        checkm2 = config["aviary_checkm2"],
        conda_prefix = config["conda_prefix"] if config["conda_prefix"] else ".",
        singlem_metapackage = config["singlem_metapackage"],
        fast = "--workflow recover_mags_no_singlem --skip-binners maxbin concoct rosella --skip-abundances --refinery-max-iterations 0" if config["aviary_speed"] == FAST_AVIARY_MODE else "",
        snakemake_profile = f"--snakemake-profile {config['snakemake_profile']}" if config["snakemake_profile"] else "",
        cluster_retries = f"--cluster-retries {config['cluster_retries']}" if config["cluster_retries"] else "",
        tmpdir = config["tmpdir"],
    localrule: True
    threads:
        int(config["aviary_threads"])//2
    resources:
        mem_mb = int(config["aviary_memory"])*1000//2,
        mem_gb = int(config["aviary_memory"])//2,
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
        "SINGLEM_METAPACKAGE_PATH={params.singlem_metapackage} "
        "TMPDIR={params.tmpdir} "
        "aviary recover "
        "--assembly {input.assembly} "
        "-1 {params.reads_1} "
        "-2 {params.reads_2} "
        "--output {params.output} "
        "{params.fast} "
        "-n {threads} "
        "-t {threads} "
        "-m {resources.mem_gb} "
        "--skip-qc "
        "{params.snakemake_profile} "
        "{params.cluster_retries} "
        "{params.dryrun} "
        "&> {log} "
        "&& touch {output} "

rule aviary_combine:
    input:
        get_coassemblies,
        elusive_clusters = output_dir + "/target/elusive_clusters.tsv",
        mapping = output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
    output:
        output_dir + "/commands/done",
    localrule: True
    shell:
        "touch {output} "
