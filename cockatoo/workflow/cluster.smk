#############
### Setup ###
#############
import os

output_dir = os.path.abspath("cluster")
logs_dir = output_dir + "/logs"

rule all:
    input:
        output_dir + "/target/elusive_edges.tsv",
        output_dir + "/target/targets.tsv",
        output_dir + "/target/elusive_clusters.tsv"

#####################
### SingleM reads ###
#####################
rule count_bp_reads:
    input:
        config["reads_1"].values()
    output:
        output_dir + "/read_size.csv"
    params:
        names=list(config["reads_1"].keys()),
        cat="zcat" if [r for r in config["reads_1"].values()][0].endswith(".gz") else "cat",
    threads:
        8
    shell:
        "parallel -k -j {threads} "
        "echo -n {{1}}, '&&' "
        "{params.cat} {{2}} '|' sed -n 2~4p '|' tr -d '\"\n\"' '|' wc -m "
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
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem pipe "
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
        singlem_metapackage=config["singlem_metapackage"]
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

######################
### Target elusive ###
######################
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
        num_coassembly_samples=config["num_coassembly_samples"],
        max_coassembly_samples=config["max_coassembly_samples"],
        max_recovery_samples=config["max_recovery_samples"],
    script:
        "scripts/cluster_graph.py"
