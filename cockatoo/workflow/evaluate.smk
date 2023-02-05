#############
### Setup ###
#############
import pandas as pd
import os

output_dir = os.path.abspath("evaluate")
logs_dir = output_dir + "/logs"

bins = config["recovered_bins"]

rule all:
    input:
        output_dir + "/evaluate/matched_hits.tsv",
        output_dir + "/evaluate/novel_hits.tsv",
        output_dir + "/evaluate/summary_stats.tsv",

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

################
### Evaluate ###
################
rule evaluate:
    input:
        recovered_otu_table = output_dir + "/summarise/bins_summarised.otu_table.tsv"
    output:
        matched_hits = output_dir + "/evaluate/matched_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
    params:
        unbinned_otu_table=config["targets"],
        binned_otu_table=config["binned"],
        elusive_edges=config["elusive_edges"],
        elusive_clusters=config["elusive_clusters"],
    script:
        "scripts/evaluate.py"

rule evaluate_plots:
    input:
        matched_hits = output_dir + "/evaluate/matched_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
    params:
        coassemble_summary=config["coassemble_summary"],
    output:
        plots_dir = directory(output_dir + "/evaluate/plots"),
        summary_stats = output_dir + "/evaluate/summary_stats.tsv",
        summary_table = output_dir + "/evaluate/summary_table.png",
    conda:
        "env/r.yml"
    script:
        "scripts/evaluate.R"
