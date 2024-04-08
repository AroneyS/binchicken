---
title: Bin chicken coassemble
---
# binchicken coassemble

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).

```bash
# Example: cluster reads into proposed coassemblies
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ...

# Example: cluster reads into proposed coassemblies based on unbinned sequences
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: cluster reads into proposed coassemblies based on unbinned sequences and coassemble only unbinned reads
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --assemble-unmapped

# Example: cluster reads into proposed coassemblies based on unbinned sequences from a specific taxa
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --taxa-of-interest "p__Planctomycetota"

# Example: find relevant samples for differential coverage binning (no coassembly)
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --single-assembly

# Example: run proposed coassemblies through aviary with cluster submission
# Create snakemake profile at ~/.config/snakemake/qsub with cluster, cluster-status, cluster-cancel, etc.
# See https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
binchicken coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --run-aviary \
  --snakemake-profile qsub --local-cores 64 --cores 64
```

Important options:

- Minimum and maximum cluster sizes can be specified (`--num-coassembly-samples` and `--max-coassembly-samples`, both default to 2)
- Maximum number of recovery samples for differential-abundance binning can be specified (`--max-recovery-samples`, default 20)
- Genomes can be provided and matching marker genes will be excluded (`--genomes`)
- Reads can be mapped to the matched bins with only unmapped reads being assembled (`--assemble-unmapped`).
- Assembly and recovery running options:
  - Run directly through Aviary (`--run-aviary`)
  - Run Aviary commands manually (see `coassemble/commands` in output)
  - Run coassemblies with differential-abudance-binning samples with the tool of your choice (see `coassemble/target/elusive_clusters.tsv` in output)
- The taxa of the considered sequences can be filtered to target a specific taxon (e.g. `--taxa-of-interest "p__Planctomycetota"`).
- Differential-abundance binning samples for single-assembly can also be found (`--single-assembly`)
- Snakemake profiles can be used to automatically submit jobs to HPC clusters (`--snakemake-profile`)

Paired end reads of form reads_1.1.fq, reads_1_1.fq and reads_1_R1.fq, where reads_1 is the sample name are automatically detected and matched to their basename.
Most intermediate files can be provided to skip intermediate steps (e.g. SingleM otu tables, read sizes or genome transcripts; see `binchicken coassemble --full-help`).

# OPTIONS

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
    and sets \--max- coassembly-size to None

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
    \--num- coassembly-samples]

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

**\--run-aviary**

  Run Aviary commands for all identified coassemblies (unless specific
    coassemblies are chosen with \--coassemblies) [default: do not]

**\--aviary-speed** {fast,comprehensive}

  Run Aviary recover in \'fast\' or \'comprehensive\' mode. Fast mode
    skips slow binners and refinement steps. [default: fast]

**\--assembly-strategy** {dynamic,metaspades,megahit}

  Assembly strategy to use with Aviary. [default: dynamic; attempts
    metaspades and if fails, switches to megahit]

**\--aviary-gtdbtk-db** *AVIARY_GTDBTK_DB*

  Path to GTDB-Tk database directory for Aviary. [default: use path
    from GTDBTK_DATA_PATH env variable]

**\--aviary-checkm2-db** *AVIARY_CHECKM2_DB*

  Path to CheckM2 database directory for Aviary. [default: use path
    from CHECKM2DB env variable]

**\--aviary-assemble-cores** *AVIARY_ASSEMBLE_CORES*

  Maximum number of cores for Aviary assemble to use. [default: 64]

**\--aviary-assemble-memory** *AVIARY_ASSEMBLE_MEMORY*

  Maximum amount of memory for Aviary assemble to use (Gigabytes).
    [default: 500]

**\--aviary-recover-cores** *AVIARY_RECOVER_CORES*

  Maximum number of cores for Aviary recover to use. [default: 32]

**\--aviary-recover-memory** *AVIARY_RECOVER_MEMORY*

  Maximum amount of memory for Aviary recover to use (Gigabytes).
    [default: 250]

# GENERAL OPTIONS

**\--output** *OUTPUT*

  Output directory [default: .]

**\--conda-prefix** *CONDA_PREFIX*

  Path to conda environment install location. [default: Use path from
    CONDA_ENV_PATH env variable]

**\--cores** *CORES*

  Maximum number of cores to use [default: 1]

**\--dryrun**

  dry run workflow

**\--snakemake-profile** *SNAKEMAKE_PROFILE*

  Snakemake profile (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cli.html#profiles).
    Can be used to submit rules as jobs to cluster engine (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cluster.html).

**\--local-cores** *LOCAL_CORES*

  Maximum number of cores to use on localrules when running in cluster
    mode [default: 1]

**\--cluster-retries** *CLUSTER_RETRIES*

  Number of times to retry a failed job when using cluster submission
    (see \`\--snakemake-profile\`) [default: 3].

**\--snakemake-args** *SNAKEMAKE_ARGS*

  Additional commands to be supplied to snakemake in the form of a
    space- prefixed single string e.g. \" \--quiet\"

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

cluster reads into proposed coassemblies

  **\$ binchicken coassemble \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \...**

cluster reads into proposed coassemblies based on unbinned sequences

  **\$ binchicken coassemble \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \... \--genomes genome_1.fna \...**

cluster reads into proposed coassemblies based on unbinned sequences and coassemble only unbinned reads

  **\$ binchicken coassemble \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \... \--genomes genome_1.fna \...
    \--assemble-unmapped**

cluster reads into proposed coassemblies based on unbinned sequences from a specific taxa

  **\$ binchicken coassemble \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \... \--genomes genome_1.fna \... \--taxa-of-interest
    \"p\_\_Planctomycetota\"**

find relevant samples for differential coverage binning (no coassembly)

  **\$ binchicken coassemble \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \... \--single-assembly**
