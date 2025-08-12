#############
### Setup ###
#############
import os
import re
from binchicken.common import pixi_run
os.umask(0o002)

output_dir = os.path.abspath("evaluate")
logs_dir = output_dir + "/logs"
scripts_dir = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')

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
        fna = output_dir + "/transcripts/{bin}_transcripts.fna",
        faa = output_dir + "/transcripts/{bin}_transcripts.faa",
    log:
        logs_dir + "/transcripts/{bin}.attempt{attempt}.log"
    params:
        prodigal_meta = "-p meta" if config["prodigal_meta"] else ""
    shell:
        f"{pixi_run} -e prodigal "
        "prodigal -i {input} -d {output.fna} -a {output.faa} "
        "{params.prodigal_meta} "
        "&> {log}"

rule singlem_pipe_bins:
    input:
        output_dir + "/transcripts/{bin}_transcripts.fna"
    output:
        output_dir + "/pipe/{bin}.otu_table.tsv"
    log:
        logs_dir + "/pipe/{bin}.attempt{attempt}.log"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    shell:
        f"{pixi_run} -e singlem "
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
        logs_dir + "/summarise/bins.attempt{attempt}.log"
    params:
        singlem_metapackage = config["singlem_metapackage"]
    shell:
        f"{pixi_run} -e singlem "
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
        logs_dir + "/cluster/original.attempt{attempt}.log"
    params:
        genomes = " ".join(config["original_bins"]),
        ani = config["cluster"],
    threads: 64
    shell:
        f"{pixi_run} -e coverm "
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
        logs_dir + "/cluster/{coassembly}.attempt{attempt}.log"
    params:
        genomes = get_cluster_genomes,
        ani = config["cluster"],
    threads: 64
    shell:
        f"{pixi_run} -e coverm "
        "coverm cluster "
        "--genome-fasta-files {params.genomes} "
        "--output-cluster-definition {output} "
        "--precluster-method finch "
        "--ani {params.ani} "
        "--threads {threads} "
        "&> {log}"

rule summarise_clusters:
    input:
        original = output_dir + "/cluster/original.txt",
        updated = expand(output_dir + "/cluster/{coassembly}.txt", coassembly=coassemblies),
    output:
        output_dir + "/evaluate/cluster_stats.csv",
    params:
        names = ["original"] + coassemblies
    shell:
        f"{pixi_run} -e general "
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
        script = scripts_dir + "/evaluate.py",
        target_otu_table = config["targets"],
        binned_otu_table = config["binned"],
        elusive_edges = config["elusive_edges"],
        elusive_clusters = config["elusive_clusters"],
        recovered_bins = config["recovered_bins"],
    threads:
        64
    log:
        logs_dir + "/evaluate/evaluate.attempt{attempt}.log"
    shell:
        f"{pixi_run} "
        "python3 {params.script} "
        "--target-otu-table {params.target_otu_table} "
        "--binned-otu-table {params.binned_otu_table} "
        "--elusive-clusters {params.elusive_clusters} "
        "--elusive-edges {params.elusive_edges} "
        "--recovered-otu-table {input.recovered_otu_table} "
        "--recovered-bins '{params.recovered_bins}' "
        "--matched-hits {output.matched_hits} "
        "--novel-hits {output.novel_hits} "
        "--summary-stats {output.summary_stats} "
        "--threads {threads} "
        "--log {log}"

rule evaluate_plots:
    input:
        matched_hits = output_dir + "/evaluate/matched_hits.tsv",
        novel_hits = output_dir + "/evaluate/novel_hits.tsv",
        cluster_summary = output_dir + "/evaluate/cluster_stats.csv" if config["cluster"] else [],
        summary_stats = output_dir + "/evaluate/summary_stats.tsv",
    params:
        script = scripts_dir + "/evaluate.R",
        cluster_summary = output_dir + "/evaluate/cluster_stats.csv" if config["cluster"] else "NULL",
        coassemble_summary = config["coassemble_summary"],
        test = config["test"] if config["test"] else "FALSE",
    output:
        plots_dir = directory(output_dir + "/evaluate/plots"),
        summary_table = output_dir + "/evaluate/summary_table.png",
    shell:
        f"{pixi_run} -e r "
        "Rscript {params.script} "
        "--matched-hits {input.matched_hits} "
        "--novel-hits {input.novel_hits} "
        "--cluster-summary {params.cluster_summary} "
        "--summary-stats {input.summary_stats} "
        "--coassemble-summary {params.coassemble_summary} "
        "--plots-dir {output.plots_dir} "
        "--summary-table {output.summary_table} "
        "--test {params.test}"
