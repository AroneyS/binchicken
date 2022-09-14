#############
### Setup ###
#############
import pandas as pd
import os

output_dir = os.path.abspath(str(config["output_subdir"]))
logs_dir = output_dir + "/logs"

if not "appraise_sequence_identity" in config: config["appraise_sequence_identity"] = 0.89
if not "min_coassembly_identity" in config:    config["min_coassembly_identity"] = 0.99
if not "taxa_of_interest" in config:           config["taxa_of_interest"] = ""
if not "min_completeness" in config:           config["min_completeness"] = 70
if not "max_contamination" in config:          config["max_contamination"] = 10

coassembly_bins = {}

if config["checkm_version"] == 1:
    completeness_col = "Completeness (CheckM1)"
    contamination_col = "Contamination (CheckM1)"
elif config["checkm_version"] == 2:
    completeness_col = "Completeness (CheckM2)"
    contamination_col = "Contamination (CheckM2)"
else:
    raise ValueError("Invalid CheckM version")

for coassembly in config["checkm_out"]:
    checkm_out = pd.read_csv(config["checkm_out"][coassembly], sep = "\t")
    passed_bins = checkm_out[(checkm_out[completeness_col] >= config["min_completeness"]) & (checkm_out[contamination_col] <= config["max_contamination"])]["Bin Id"].to_list()
    coassembly_bins[coassembly] = passed_bins

bins = {"-".join([c, b]): os.path.join(config["recovered_bins"][c], b + ".fna") for c in coassembly_bins for b in coassembly_bins[c]}

rule all:
    input:
        output_dir + "/evaluate/unbinned_hits.tsv",
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
        "/home/aroneys/conda_env_yamls/prodigal.yaml"
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
        "/home/aroneys/src/singlem/singlem.yml"
    shell:
        "/home/aroneys/bin/singlem-dev pipe "
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
        "/home/aroneys/src/singlem/singlem.yml"
    shell:
        "/home/aroneys/bin/singlem-dev summarise "
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
        unbinned_hits = output_dir + "/evaluate/unbinned_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
    params:
        unbinned_otu_table=config["targets"],
        elusive_edges=config["elusive_edges"],
        elusive_clusters=config["elusive_clusters"],
    script:
        "scripts/evaluate.py"

rule evaluate_plots:
    input:
        unbinned_hits = output_dir + "/evaluate/unbinned_hits.tsv",
    output:
        plots_dir = directory(output_dir + "/evaluate/plots"),
        summary_stats = output_dir + "/evaluate/summary_stats.tsv",
    conda:
        "env/r.yml"
    script:
        "scripts/evaluate.R"

# Join new hits with old target names by sequence
# Check for novelty (if too many then cluster)
# Check for expected/unexpected elusive targets
# Summarise taxonomy of total recovered vs newly binned (expected/unexpected)
