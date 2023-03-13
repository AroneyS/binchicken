#############
### Setup ###
#############
import pandas as pd
import os
import re

output_dir = os.path.abspath("evaluate")
logs_dir = output_dir + "/logs"

bins = config["recovered_bins"]
coassemblies = list(set([re.search(r"coassembly_\d+", b)[0] for b in bins.keys()]))

def get_cluster_genomes(wildcards):
    coassembly = wildcards.coassembly
    new_bins = [bins[b] for b in bins if coassembly in b]

    return " ".join(config["original_bins"] + new_bins)

rule all:
    input:
        output_dir + "/evaluate/matched_hits.tsv",
        output_dir + "/evaluate/novel_hits.tsv",
        output_dir + "/evaluate/summary_stats.tsv",
        output_dir + "/evaluate/summary_table.png",

######################
### Recovered bins ###
######################
rule prodigal_bins:
    input:
        lambda wildcards: bins[wildcards.bin]
    output:
        output_dir + "/transcripts/{bin}_transcripts.fna"
    log:
        logs_dir + "/transcripts/{bin}.log"
    conda:
        "env/prodigal.yml"
    shell:
        "prodigal -i {input} -d {output} "
        "&> {log}"

rule singlem_pipe_bins:
    input:
        output_dir + "/transcripts/{bin}_transcripts.fna"
    output:
        output_dir + "/pipe/{bin}.otu_table.tsv"
    log:
        logs_dir + "/pipe/{bin}.log"
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
        expand(output_dir + "/pipe/{bin}.otu_table.tsv", bin=bins.keys())
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

###############
### Cluster ###
###############
rule cluster_original_bins:
    output:
        output_dir + "/cluster/original.txt",
    log:
        logs_dir + "/cluster/original.log"
    params:
        genomes=" ".join(config["original_bins"]),
        ani=config["cluster"],
    threads: 64
    conda:
        "env/coverm.yml"
    shell:
        "coverm cluster "
        "--genome-fasta-files {params.genomes} "
        "--output-cluster-definition {output} "
        "--precluster-method finch "
        "--ani {params.ani} "
        "--threads {threads} "
        "&> {log}"

rule cluster_updated_bins:
    output:
        output_dir + "/cluster/{coassembly}.txt",
    log:
        logs_dir + "/cluster/{coassembly}.log"
    params:
        genomes=get_cluster_genomes,
        ani=config["cluster"],
    threads: 64
    conda:
        "env/coverm.yml"
    shell:
        "coverm cluster "
        "--genome-fasta-files {params.genomes} "
        "--output-cluster-definition {output} "
        "--precluster-method finch "
        "--ani {params.ani} "
        "--threads {threads} "
        "&> {log}"

rule summarise_clusters:
    input:
        original=output_dir + "/cluster/original.txt",
        updated=expand(output_dir + "/cluster/{coassembly}.txt", coassembly=coassemblies),
    output:
        output_dir + "/evaluate/cluster_stats.csv",
    params:
        names=["original"] + coassemblies
    shell:
        "parallel -k "
        "echo -n {{1}}, '&&' "
        "cut -f1 {{2}} '|' uniq '|' wc -l "
        "::: {params.names} :::+ {input.original} {input.updated} "
        "> {output}"

################
### Evaluate ###
################
rule evaluate:
    input:
        recovered_otu_table = output_dir + "/summarise/bins_summarised.otu_table.tsv"
    output:
        matched_hits = output_dir + "/evaluate/matched_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
        summary_stats = output_dir + "/evaluate/summary_stats.tsv",
    params:
        unbinned_otu_table=config["targets"],
        binned_otu_table=config["binned"],
        elusive_edges=config["elusive_edges"],
        elusive_clusters=config["elusive_clusters"],
        recovered_bins=config["recovered_bins"],
    threads:
        64
    script:
        "scripts/evaluate.py"

rule evaluate_plots:
    input:
        matched_hits = output_dir + "/evaluate/matched_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
        cluster_summary = output_dir + "/evaluate/cluster_stats.csv" if config["cluster"] else [],
        summary_stats = output_dir + "/evaluate/summary_stats.tsv",
    params:
        coassemble_summary=config["coassemble_summary"],
    output:
        plots_dir = directory(output_dir + "/evaluate/plots"),
        summary_table = output_dir + "/evaluate/summary_table.png",
    conda:
        "env/r.yml"
    script:
        "scripts/evaluate.R"
