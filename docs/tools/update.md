---
title: Bin chicken update
---
# binchicken update

Applies further processing to a previous Bin chicken coassemble run.

Any combinations of the following:

- Generating unmapped reads files (`--assemble-unmapped`)
- Running assembly/recovery through Aviary (`--run-aviary`)
- Downloading SRA reads (`--sra`)

```bash
# Example: update previous run to perform unmapping
binchicken update --coassemble-output coassemble_dir --assemble-unmapped --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to run specific coassemblies
binchicken update --coassemble-output coassemble_dir --run-aviary --coassemblies coassembly_0 ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to download SRA reads
binchicken update --coassemble-output coassemble_dir --sra --forward SRA000001 ... --genomes genome_1.fna ...
```

# OPTIONS

# INPUT ARGUMENTS

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

**\--sra**

  Download reads from SRA (read argument still required). Also sets
    \--run-qc.

# COASSEMBLY OPTIONS

**\--coassemble-output** *COASSEMBLE_OUTPUT*

  Output dir from coassemble subcommand

**\--coassemble-unbinned** *COASSEMBLE_UNBINNED*

  SingleM appraise unbinned output from Bin chicken coassemble
    (alternative to \--coassemble-output)

**\--coassemble-binned** *COASSEMBLE_BINNED*

  SingleM appraise binned output from Bin chicken coassemble
    (alternative to \--coassemble-output)

**\--coassemble-targets** *COASSEMBLE_TARGETS*

  Target sequences output from Bin chicken coassemble (alternative to
    \--coassemble-output)

**\--coassemble-elusive-edges** *COASSEMBLE_ELUSIVE_EDGES*

  Elusive edges output from Bin chicken coassemble (alternative to
    \--coassemble-output)

**\--coassemble-elusive-clusters** *COASSEMBLE_ELUSIVE_CLUSTERS*

  Elusive clusters output from Bin chicken coassemble (alternative to
    \--coassemble-output)

**\--coassemble-summary** *COASSEMBLE_SUMMARY*

  Summary output from Bin chicken coassemble (alternative to
    \--coassemble-output)

**\--coassemblies** *COASSEMBLIES* [*COASSEMBLIES* \...]

  Choose specific coassemblies from elusive clusters (e.g.
    coassembly_0)

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

update previous run to perform unmapping

  **\$ binchicken update \--coassemble-output coassemble_dir
    \--assemble-unmapped \--forward reads_1.1.fq \... \--reverse
    reads_1.2.fq \... \--genomes genome_1.fna \...**

update previous run to run specific coassemblies

  **\$ binchicken update \--coassemble-output coassemble_dir
    \--run-aviary \--coassemblies coassembly_0 \... \--forward
    reads_1.1.fq \... \--reverse reads_1.2.fq \... \--genomes
    genome_1.fna \...**

update previous run to download SRA reads

  **\$ binchicken update \--coassemble-output coassemble_dir \--sra
    \--forward SRA000001 \... \--genomes genome_1.fna \...**
