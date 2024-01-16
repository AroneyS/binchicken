---
title: Bin chicken evaluate
---
# binchicken evaluate

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

```bash
# Example: evaluate a completed coassembly
binchicken evaluate --coassemble-output coassemble_dir --aviary-outputs coassembly_0_dir ...

# Example: evaluate a completed coassembly by providing genomes directly
binchicken evaluate --coassemble-output coassemble_dir --new-genomes genome_1.fna ... --coassembly-run coassembly_0
```

Defaults to using genomes (from the provided coassemble outputs) with at least 70% complete and at most 10% contamination as estimated by CheckM2.

# OPTIONS

# BASE INPUT ARGUMENTS

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

**\--aviary-outputs** *AVIARY_OUTPUTS* [*AVIARY_OUTPUTS* \...]

  Output dir from Aviary coassembly and recover commands produced by
    coassemble subcommand

**\--new-genomes** *NEW_GENOMES* [*NEW_GENOMES* \...]

  New genomes to evaluate (alternative to \--aviary-outputs, also
    requires \--coassembly-run)

**\--new-genomes-list** *NEW_GENOMES_LIST*

  New genomes to evaluate (alternative to \--aviary-outputs, also
    requires \--coassembly-run) newline separated

**\--coassembly-run** *COASSEMBLY_RUN*

  Name of coassembly run to produce new genomes (alternative to
    \--aviary-outputs, also requires \--new-genomes)

**\--singlem-metapackage** *SINGLEM_METAPACKAGE*

  SingleM metapackage for sequence searching

**\--prodigal-meta**

  Use prodigal \"-p meta\" argument (for testing)

# EVALUATION OPTIONS

**\--checkm-version** *CHECKM_VERSION*

  CheckM version to use to quality cutoffs [default: 2]

**\--min-completeness** *MIN_COMPLETENESS*

  Include bins with at least this minimum completeness [default: 70]

**\--max-contamination** *MAX_CONTAMINATION*

  Include bins with at most this maximum contamination [default: 10]

# CLUSTER OPTIONS

**\--cluster**

  Cluster new and original genomes and report number of new clusters

**\--cluster-ani** *CLUSTER_ANI*

  Cluster using this sequence identity [default: 86%]

**\--genomes** *GENOMES* [*GENOMES* \...]

  Original genomes used as references for coassemble subcommand

**\--genomes-list** *GENOMES_LIST*

  Original genomes used as references for coassemble subcommand
    newline separated

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

evaluate a completed coassembly

  **\$ binchicken evaluate \--coassemble-output coassemble_dir
    \--aviary-outputs coassembly_0_dir \...**

evaluate a completed coassembly by providing genomes directly

  **\$ binchicken evaluate \--coassemble-output coassemble_dir
    \--new-genomes genome_1.fna \... \--coassembly-run coassembly_0**
