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
# Example: update previous run to run specific coassemblies
binchicken update --coassemble-output coassemble_dir --run-aviary \
    --coassemblies coassembly_0 ...

# Example: update previous run to perform unmapping
binchicken update --coassemble-output coassemble_dir --assemble-unmapped

# Example: update previous run to download SRA reads
# Note: requires sample names to be SRA IDs (e.g. SRA123456)
binchicken update --coassemble-output coassemble_dir --sra

# Example: update previous run to download SRA reads, perform unmapping and run specific coassemblies
binchicken update --coassemble-output coassemble_dir --sra \
    --assemble-unmapped \
    --run-aviary --coassemblies coassembly_0 ...
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
    \--coassemble- output)

**\--coassemble-elusive-clusters** *COASSEMBLE_ELUSIVE_CLUSTERS*

  Elusive clusters output from Bin chicken coassemble (alternative to
    \--coassemble-output)

**\--coassemble-summary** *COASSEMBLE_SUMMARY*

  Summary output from Bin chicken coassemble (alternative to
    \--coassemble- output)

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

update previous run to run specific coassemblies

  **\$ binchicken update \--coassemble-output coassemble_dir
    \--run-aviary \--coassemblies coassembly_0 \...**

update previous run to perform unmapping

  **\$ binchicken update \--coassemble-output coassemble_dir
    \--assemble-unmapped**

update previous run to download SRA reads

  **\$ binchicken update \--coassemble-output coassemble_dir \--sra**

update previous run to download SRA reads, perform unmapping and run specific coassemblies

  **\$ binchicken update \--coassemble-output coassemble_dir \--sra
    \--assemble-unmapped \--run-aviary \--coassemblies coassembly_0
    \...**
