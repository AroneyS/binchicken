#############
### Setup ###
#############
ruleorder: singlem_appraise_gzip_archive > singlem_appraise

import os

output_dir = os.path.abspath(str(config["output_subdir"]))
logs_dir = output_dir + "/logs"

singlem_bin = "/home/aroneys/bin/singlem-dev"
singlem_conda = "/home/aroneys/src/singlem/singlem.yml"

if not "min_contig_size" in config:            config["min_contig_size"] = 2500
if not "appraise_sequence_identity" in config: config["appraise_sequence_identity"] = 0.89
if not "max_coassembly_size" in config:        config["max_coassembly_size"] = 50
if not "max_coassembly_samples" in config:     config["max_coassembly_samples"] = 5
if not "min_coassembly_identity" in config:    config["min_coassembly_identity"] = 0.99
if not "min_coassembly_coverage" in config:    config["min_coassembly_coverage"] = 10
if not "max_recovery_samples" in config:       config["max_recovery_samples"] = 20
if not "taxa_of_interest" in config:           config["taxa_of_interest"] = ""

rule all:
    input:
        output_dir + "/target/coassemble_commands.sh",
        output_dir + "/target/recover_commands.sh"

#####################
### SingleM reads ###
#####################
rule count_bp_reads:
    input:
        config["reads_1"].values()
    output:
        output_dir + "/read_size.csv"
    params:
        names=list(config["reads_1"].keys())
    threads:
        8
    shell:
        "parallel -k -j {threads} "
        "echo -n {{1}}, '&&' "
        "zcat {{2}} '|' sed -n 2~4p '|' tr -d '\"\n\"' '|' wc -m "
        "::: {params.names} :::+ {input} "
        "> {output}"

rule singlem_pipe_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read],
    output:
        output_dir + "/pipe/{read}_read.otu_table.tsv"
    log:
        logs_dir + "/pipe/{read}_read.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} pipe "
        "--forward {input.reads_1} "
        "--reverse {input.reads_2} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

rule singlem_summarise_reads:
    input:
        output_dir + "/pipe/{read}_read.otu_table.tsv"
    output:
        output_dir + "/summarise/{read}_summarised.otu_table.tsv"
    log:
        logs_dir + "/summarise/{read}.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--exclude-off-target-hits "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

####################
### SingleM bins ###
####################
rule singlem_pipe_bins:
    input:
        lambda wildcards: config["bin_transcripts"][wildcards.bin]
    output:
        output_dir + "/pipe/{bin}_bin.otu_table.tsv"
    log:
        logs_dir + "/pipe/{bin}_bin.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} pipe "
        "--forward {input} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

rule singlem_summarise_bins:
    input:
        expand(output_dir + "/pipe/{bin}_bin.otu_table.tsv", bin=config["bin_transcripts"])
    output:
        output_dir + "/summarise/bins_summarised.otu_table.tsv"
    log:
        logs_dir + "/summarise/bins.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--exclude-off-target-hits "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

##################################
### SingleM appraise/summarise ###
##################################
rule singlem_appraise:
    input:
        reads=output_dir + "/summarise/{read}_summarised.otu_table.tsv",
        bins=output_dir + "/summarise/bins_summarised.otu_table.tsv",
    output:
        unbinned=output_dir + "/appraise/{read}_unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/{read}_binned.otu_table.tsv",
    log:
        logs_dir + "/appraise/{read}.log"
    params:
        singlem_bin = singlem_bin,
        sequence_identity=config["appraise_sequence_identity"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} appraise "
        "--metagenome-otu-tables {input.reads} "
        "--genome-otu-tables {input.bins} "
        "--output-unaccounted-for-otu-table {output.unbinned} "
        "--output-binned-otu-table {output.binned} "
        "--imperfect "
        "--sequence-identity {params.sequence_identity} "
        "&> {log}"

rule singlem_appraise_gzip_archive:
    input:
        reads=lambda wildcards: config["singlem_gzip_archive"][wildcards.read],
        bins=output_dir + "/summarise/bins_summarised.otu_table.tsv",
    output:
        unbinned=output_dir + "/appraise/{read}_unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/{read}_binned.otu_table.tsv",
    log:
        logs_dir + "/appraise/{read}.log"
    params:
        singlem_bin = singlem_bin,
        sequence_identity=config["appraise_sequence_identity"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} appraise "
        "--metagenome-gzip-archive-otu-tables {input.reads} "
        "--genome-otu-tables {input.bins} "
        "--output-unaccounted-for-otu-table {output.unbinned} "
        "--output-binned-otu-table {output.binned} "
        "--imperfect "
        "--sequence-identity {params.sequence_identity} "
        "&> {log}"

rule singlem_summarise_unbinned:
    input:
        expand(output_dir + "/appraise/{read}_unbinned.otu_table.tsv", read=config["reads_1"])
    output:
        output_dir + "/summarise/unbinned.otu_table.tsv"
    log:
        logs_dir + "/summarise/unbinned.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

rule singlem_summarise_binned:
    input:
        expand(output_dir + "/appraise/{read}_binned.otu_table.tsv", read=config["reads_1"])
    output:
        output_dir + "/summarise/binned.otu_table.tsv"
    log:
        logs_dir + "/summarise/binned.log"
    params:
        singlem_bin = singlem_bin,
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        singlem_conda
    shell:
        "{params.singlem_bin} summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

######################
### Target elusive ###
######################
rule target_elusive:
    input:
        unbinned=output_dir + "/summarise/unbinned.otu_table.tsv"
    output:
        output_edges=output_dir + "/target/elusive_edges.tsv",
        output_targets=output_dir + "/target/targets.tsv",
    log:
        logs_dir + "/elusive_targets.log"
    params:
        min_coassembly_coverage=config["min_coassembly_coverage"],
        taxa_of_interest=config["taxa_of_interest"],
    script:
        "scripts/target_elusive.py"

rule cluster_graph:
    input:
        elusive_edges=output_dir + "/target/elusive_edges.tsv",
        read_size=output_dir + "/read_size.csv",
    output:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv"
    log:
        logs_dir + "/cluster_graph.log"
    params:
        max_coassembly_size=config["max_coassembly_size"],
        max_coassembly_samples=config["max_coassembly_samples"],
        max_recovery_samples=config["max_recovery_samples"],
    script:
        "scripts/cluster_graph.py"

rule coassemble_commands:
    input:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv"
    output:
        coassemble_commands=output_dir + "/target/coassemble_commands.sh",
        recover_commands=output_dir + "/target/recover_commands.sh"
    log:
        logs_dir + "/coassemble_commands.log"
    script:
        "scripts/coassemble_commands.py"

# Co-assemble identified samples
# Provide samples together space-separated to aviary (will automatically use megahit)
