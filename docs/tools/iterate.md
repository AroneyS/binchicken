---
title: Bin chicken iterate
---
# binchicken iterate

Run a further iteration of coassemble, including newly recovered bins.
All [coassemble](/tools/coassemble) options are available and can be altered for the new iteration (see [example workflow](/)).

```bash
# Example: rerun coassemble, adding new bins to database
binchicken iterate --coassemble-output coassemble_dir

# Example: rerun coassemble, adding new bins to database, providing genomes directly
binchicken iterate --coassemble-output coassemble_dir --new-genomes new_genome_1.fna
```

Defaults to using genomes (from the provided coassemble outputs) with at least 70% complete and at most 10% contamination as estimated by CheckM2.
Automatically excludes previous coassemblies.

# OPTIONS

# ITERATION OPTIONS

**\--iteration** *ITERATION*

  Iteration number used for unique bin naming

**\--aviary-outputs** *AVIARY_OUTPUTS* [*AVIARY_OUTPUTS* \...]

  Output dir from Aviary coassembly and recover commands produced by
    coassemble subcommand

**\--new-genomes** *NEW_GENOMES* [*NEW_GENOMES* \...]

  New genomes to iterate (alternative to \--aviary-outputs)

**\--new-genomes-list** *NEW_GENOMES_LIST*

  New genomes to iterate (alternative to \--aviary-outputs) newline
    separated

**\--new-genome-singlem** *NEW_GENOME_SINGLEM*

  Combined SingleM otu tables for new genome transcripts. If provided,
    genome SingleM is skipped

**\--elusive-clusters** *ELUSIVE_CLUSTERS* [*ELUSIVE_CLUSTERS* \...]

  Previous elusive_clusters.tsv files produced by coassemble
    subcommand (used to check for duplicated coassembly suggestions)

**\--coassemble-output** *COASSEMBLE_OUTPUT*

  Output dir from coassemble subcommand

**\--coassemble-unbinned** *COASSEMBLE_UNBINNED*

  SingleM appraise unbinned output from Bin chicken coassemble
    (alternative to \--coassemble-output)

**\--coassemble-binned** *COASSEMBLE_BINNED*

  SingleM appraise binned output from Bin chicken coassemble
    (alternative to \--coassemble-output)

**\--checkm-version** *CHECKM_VERSION*

  CheckM version to use to quality cutoffs [default: 2]

**\--min-completeness** *MIN_COMPLETENESS*

  Include bins with at least this minimum completeness [default: 70]

**\--max-contamination** *MAX_CONTAMINATION*

  Include bins with at most this maximum contamination [default: 10]

# BASE INPUT ARGUMENTS

**\--forward**, **\--reads**, **\--sequences** *FORWARD* [*FORWARD* \...]

  input forward/unpaired nucleotide read sequence(s)

**\--forward-list**, **\--reads-list**, **\--sequences-list** *FORWARD_LIST*

  input forward/unpaired nucleotide read sequence(s) newline separated

**\--reverse** *REVERSE* [*REVERSE* \...]

  input reverse nucleotide read sequence(s)

**\--reverse-list** *REVERSE_LIST*

  input reverse nucleotide read sequence(s) newline separated

**\--genomes** *GENOMES* [*GENOMES* \...]

  Reference genomes for read mapping

**\--genomes-list** *GENOMES_LIST*

  Reference genomes for read mapping newline separated

**\--coassembly-samples** *COASSEMBLY_SAMPLES* [*COASSEMBLY_SAMPLES* \...]

  Restrict coassembly to these samples. Remaining samples will still
    be used for recovery [default: use all samples]

**\--coassembly-samples-list** *COASSEMBLY_SAMPLES_LIST*

  Restrict coassembly to these samples, newline separated. Remaining
    samples will still be used for recovery [default: use all samples]

**\--singlem-metapackage** *SINGLEM_METAPACKAGE*

  SingleM metapackage for sequence searching. [default: use path from
    SINGLEM_METAPACKAGE_PATH env variable]

# INTERMEDIATE RESULTS INPUT ARGUMENTS

**\--sample-singlem** *SAMPLE_SINGLEM* [*SAMPLE_SINGLEM* \...]

  SingleM otu tables for each sample, in the form \"[sample
    name]\_read.otu_table.tsv\". If provided, SingleM pipe sample is
    skipped

**\--sample-singlem-list** *SAMPLE_SINGLEM_LIST*

  SingleM otu tables for each sample, in the form \"[sample
    name]\_read.otu_table.tsv\" newline separated. If provided, SingleM
    pipe sample is skipped

**\--sample-singlem-dir** *SAMPLE_SINGLEM_DIR*

  Directory containing SingleM otu tables for each sample, in the form
    \"[sample name]\_read.otu_table.tsv\". If provided, SingleM pipe
    sample is skipped

**\--sample-query** *SAMPLE_QUERY* [*SAMPLE_QUERY* \...]

  Queried SingleM otu tables for each sample against genome database,
    in the form \"[sample name]\_query.otu_table.tsv\". If provided,
    SingleM pipe and appraise are skipped

**\--sample-query-list** *SAMPLE_QUERY_LIST*

  Queried SingleM otu tables for each sample against genome database,
    in the form \"[sample name]\_query.otu_table.tsv\" newline
    separated. If provided, SingleM pipe and appraise are skipped

**\--sample-query-dir** *SAMPLE_QUERY_DIR*

  Directory containing Queried SingleM otu tables for each sample
    against genome database, in the form \"[sample
    name]\_query.otu_table.tsv\". If provided, SingleM pipe and
    appraise are skipped

**\--sample-read-size** *SAMPLE_READ_SIZE*

  Comma separated list of sample name and size (bp). If provided,
    sample read counting is skipped

**\--genome-transcripts** *GENOME_TRANSCRIPTS* [*GENOME_TRANSCRIPTS* \...]

  Genome transcripts for reference database, in the form
    \"[genome]\_protein.fna\"

**\--genome-transcripts-list** *GENOME_TRANSCRIPTS_LIST*

  Genome transcripts for reference database, in the form
    \"[genome]\_protein.fna\" newline separated

**\--genome-singlem** *GENOME_SINGLEM*

  Combined SingleM otu tables for genome transcripts. If provided,
    genome SingleM is skipped

# CLUSTERING OPTIONS

**\--taxa-of-interest** *TAXA_OF_INTEREST*

  Only consider sequences from this GTDB taxa (e.g.
    p\_\_Planctomycetota) [default: all]

**\--appraise-sequence-identity** *APPRAISE_SEQUENCE_IDENTITY*

  Minimum sequence identity for SingleM appraise against reference
    database [default: 86%, Genus-level]

**\--min-sequence-coverage** *MIN_SEQUENCE_COVERAGE*

  Minimum combined coverage for sequence inclusion [default: 10]

**\--single-assembly**

  Skip appraise to discover samples to differential abundance binning.
    Forces \--num-coassembly-samples and \--max-coassembly-samples to 1
    and sets \--max-coassembly-size to None

**\--exclude-coassemblies** *EXCLUDE_COASSEMBLIES* [*EXCLUDE_COASSEMBLIES* \...]

  List of coassemblies to exclude, space separated, in the form
    \"sample_1,sample_2\"

**\--exclude-coassemblies-list** *EXCLUDE_COASSEMBLIES_LIST*

  List of coassemblies to exclude, space separated, in the form
    \"sample_1,sample_2\", newline separated

**\--num-coassembly-samples** *NUM_COASSEMBLY_SAMPLES*

  Number of samples per coassembly cluster [default: 2]

**\--max-coassembly-samples** *MAX_COASSEMBLY_SAMPLES*

  Upper bound for number of samples per coassembly cluster [default:
    \--num-coassembly-samples]

**\--max-coassembly-size** *MAX_COASSEMBLY_SIZE*

  Maximum size (Gbp) of coassembly cluster [default: 50Gbp]

**\--max-recovery-samples** *MAX_RECOVERY_SAMPLES*

  Upper bound for number of related samples to use for differential
    abundance binning [default: 20]

**\--prodigal-meta**

  Use prodigal \"-p meta\" argument (for testing)

# COASSEMBLY OPTIONS

**\--assemble-unmapped**

  Only assemble reads that do not map to reference genomes

**\--run-qc**

  Run Fastp QC on reads

**\--unmapping-min-appraised** *UNMAPPING_MIN_APPRAISED*

  Minimum fraction of sequences binned to justify unmapping [default:
    0.1]

**\--unmapping-max-identity** *UNMAPPING_MAX_IDENTITY*

  Maximum sequence identity of mapped sequences kept for coassembly
    [default: 99%]

**\--unmapping-max-alignment** *UNMAPPING_MAX_ALIGNMENT*

  Maximum percent alignment of mapped sequences kept for coassembly
    [default: 99%]

**\--aviary-speed** {fast,comprehensive}

  Run Aviary recover in \'fast\' or \'comprehensive\' mode. Fast mode
    skips slow binners and refinement steps.

**\--run-aviary**

  Run Aviary commands for all identified coassemblies (unless
    specified)

**\--aviary-gtdbtk-db** *AVIARY_GTDBTK_DB*

  Path to GTDB-Tk database directory for Aviary. [default: use path
    from GTDBTK_DATA_PATH env variable]

**\--aviary-checkm2-db** *AVIARY_CHECKM2_DB*

  Path to CheckM2 database directory for Aviary. [default: use path
    from CHECKM2DB env variable]

**\--aviary-cores** *AVIARY_CORES*

  Maximum number of cores for Aviary to use. Half used for recovery.

**\--aviary-memory** *AVIARY_MEMORY*

  Maximum amount of memory for Aviary to use (Gigabytes). Half used
    for recovery

# GENERAL OPTIONS

**\--output** *OUTPUT*

  Output directory [default: .]

**\--conda-prefix** *CONDA_PREFIX*

  Path to conda environment install location. [default: Use path from
    CONDA_ENV_PATH env variable]

**\--cores** *CORES*

  Maximum number of cores to use

**\--dryrun**

  dry run workflow

**\--snakemake-profile** *SNAKEMAKE_PROFILE*

  Snakemake profile (see
    https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
    Can be used to submit rules as jobs to cluster engine (see
    https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

**\--local-cores** *LOCAL_CORES*

  Maximum number of cores to use on localrules when running in cluster
    mode

**\--cluster-retries** *CLUSTER_RETRIES*

  Number of times to retry a failed job when using cluster submission
    (see \`\--snakemake-profile\`).

**\--snakemake-args** *SNAKEMAKE_ARGS*

  Additional commands to be supplied to snakemake in the form of a
    space-prefixed single string e.g. \" \--quiet\"

**\--tmp-dir** *TMP_DIR*

  Path to temporary directory. [default: Use path from TMPDIR env
    variable]

# OTHER GENERAL OPTIONS

**\--debug**

  output debug information

**\--version**

  output version information and quit

**\--quiet**

  only output errors

**\--full-help**

  print longer help message

**\--full-help-roff**

  print longer help message in ROFF (manpage) format

# EXAMPLES

rerun coassemble, adding new bins to database

  **\$ binchicken iterate \--coassemble-output coassemble_dir**

rerun coassemble, adding new bins to database, providing genomes directly

  **\$ binchicken iterate \--coassemble-output coassemble_dir
    \--new-genomes new_genome_1.fna**
