---
title: Bin Chicken iterate
---
# binchicken iterate

Run a further iteration of coassemble, including newly recovered bins.
All [coassemble](/tools/coassemble) options are available and can be altered for the new iteration (see [example workflow](/)).

```bash
# Example: rerun coassemble, adding new bins to database
binchicken iterate --coassemble-output coassemble_dir

# Example: rerun coassemble, adding new bins to database, providing genomes directly
binchicken iterate --coassemble-output coassemble_dir \
    --new-genomes new_genome_1.fna
```

Defaults to using genomes (from the provided coassemble outputs) with at least 70% complete and at most 10% contamination as estimated by CheckM2.
Alternatively, selected genomes can be provided directly with `--new-genomes`.
Automatically excludes previous coassemblies.

# OPTIONS

# ITERATION OPTIONS

**\--iteration** *ITERATION*

  Iteration number used for unique bin naming [default: 0]

<!-- -->

**\--aviary-outputs** *AVIARY_OUTPUTS* [*AVIARY_OUTPUTS* \...]

  Output dir from Aviary coassembly and recover commands produced by
    coassemble subcommand

<!-- -->

**\--new-genomes** *NEW_GENOMES* [*NEW_GENOMES* \...]

  New genomes to iterate (alternative to \--aviary-outputs)

<!-- -->

**\--new-genomes-list** *NEW_GENOMES_LIST*

  New genomes to iterate (alternative to \--aviary-outputs) newline
    separated

<!-- -->

**\--new-genome-singlem** *NEW_GENOME_SINGLEM*

  Combined SingleM otu tables for new genome transcripts. If provided,
    genome SingleM is skipped

<!-- -->

**\--elusive-clusters** *ELUSIVE_CLUSTERS* [*ELUSIVE_CLUSTERS* \...]

  Previous elusive_clusters.tsv files produced by coassemble
    subcommand (used to check for duplicated coassembly suggestions)

<!-- -->

**\--coassemble-output** *COASSEMBLE_OUTPUT*

  Output dir from coassemble subcommand

<!-- -->

**\--coassemble-unbinned** *COASSEMBLE_UNBINNED*

  SingleM appraise unbinned output from Bin Chicken coassemble
    (alternative to \--coassemble-output)

<!-- -->

**\--coassemble-binned** *COASSEMBLE_BINNED*

  SingleM appraise binned output from Bin Chicken coassemble
    (alternative to \--coassemble-output)

<!-- -->

**\--checkm-version** *CHECKM_VERSION*

  CheckM version to use to quality cutoffs [default: 2]

<!-- -->

**\--min-completeness** *MIN_COMPLETENESS*

  Include bins with at least this minimum completeness [default: 70]

<!-- -->

**\--max-contamination** *MAX_CONTAMINATION*

  Include bins with at most this maximum contamination [default: 10]

# BASE INPUT ARGUMENTS

**\--forward**, **\--reads**, **\--sequences** *FORWARD* [*FORWARD* \...]

  input forward/unpaired nucleotide read sequence(s)

<!-- -->

**\--forward-list**, **\--reads-list**, **\--sequences-list** *FORWARD_LIST*

  input forward/unpaired nucleotide read sequence(s) newline separated

<!-- -->

**\--reverse** *REVERSE* [*REVERSE* \...]

  input reverse nucleotide read sequence(s)

<!-- -->

**\--reverse-list** *REVERSE_LIST*

  input reverse nucleotide read sequence(s) newline separated

<!-- -->

**\--genomes** *GENOMES* [*GENOMES* \...]

  Reference genomes for read mapping

<!-- -->

**\--genomes-list** *GENOMES_LIST*

  Reference genomes for read mapping newline separated

<!-- -->

**\--coassembly-samples** *COASSEMBLY_SAMPLES* [*COASSEMBLY_SAMPLES* \...]

  Restrict coassembly to these samples. Remaining samples will still
    be used for recovery [default: use all samples]

<!-- -->

**\--coassembly-samples-list** *COASSEMBLY_SAMPLES_LIST*

  Restrict coassembly to these samples, newline separated. Remaining
    samples will still be used for recovery [default: use all samples]

<!-- -->

**\--anchor-samples** *ANCHOR_SAMPLES* [*ANCHOR_SAMPLES* \...]

  Samples to use as anchors for coassembly, all coassemblies will
    contain at least one anchor sample. [default: no restriction]

<!-- -->

**\--anchor-samples-list** *ANCHOR_SAMPLES_LIST*

  Samples to use as anchors for coassembly, all coassemblies will
    contain at least one anchor sample, newline separated. [default: no
    restriction]

<!-- -->

**\--singlem-metapackage** *SINGLEM_METAPACKAGE*

  SingleM metapackage for sequence searching. [default: use path from
    SINGLEM_METAPACKAGE_PATH env variable]

# INTERMEDIATE RESULTS INPUT ARGUMENTS

**\--sample-singlem** *SAMPLE_SINGLEM* [*SAMPLE_SINGLEM* \...]

  SingleM otu tables for each sample, in the form \"[sample
    name]\_read.otu_table.tsv\". If provided, SingleM pipe sample is
    skipped

<!-- -->

**\--sample-singlem-list** *SAMPLE_SINGLEM_LIST*

  SingleM otu tables for each sample, in the form \"[sample
    name]\_read.otu_table.tsv\" newline separated. If provided, SingleM
    pipe sample is skipped

<!-- -->

**\--sample-singlem-dir** *SAMPLE_SINGLEM_DIR*

  Directory containing SingleM otu tables for each sample, in the form
    \"[sample name]\_read.otu_table.tsv\". If provided, SingleM pipe
    sample is skipped

<!-- -->

**\--sample-query** *SAMPLE_QUERY* [*SAMPLE_QUERY* \...]

  Queried SingleM otu tables for each sample against genome database,
    in the form \"[sample name]\_query.otu_table.tsv\". If provided,
    SingleM pipe and appraise are skipped

<!-- -->

**\--sample-query-list** *SAMPLE_QUERY_LIST*

  Queried SingleM otu tables for each sample against genome database,
    in the form \"[sample name]\_query.otu_table.tsv\" newline
    separated. If provided, SingleM pipe and appraise are skipped

<!-- -->

**\--sample-query-dir** *SAMPLE_QUERY_DIR*

  Directory containing Queried SingleM otu tables for each sample
    against genome database, in the form \"[sample
    name]\_query.otu_table.tsv\". If provided, SingleM pipe and
    appraise are skipped

<!-- -->

**\--sample-read-size** *SAMPLE_READ_SIZE*

  Comma separated list of sample name and size (bp). If provided,
    sample read counting is skipped

<!-- -->

**\--genome-transcripts** *GENOME_TRANSCRIPTS* [*GENOME_TRANSCRIPTS* \...]

  Genome transcripts for reference database, in the form
    \"[genome]\_protein.fna\"

<!-- -->

**\--genome-transcripts-list** *GENOME_TRANSCRIPTS_LIST*

  Genome transcripts for reference database, in the form
    \"[genome]\_protein.fna\" newline separated

<!-- -->

**\--genome-singlem** *GENOME_SINGLEM*

  Combined SingleM otu tables for genome transcripts. If provided,
    genome SingleM is skipped

# CLUSTERING OPTIONS

**\--taxa-of-interest** *TAXA_OF_INTEREST*

  Only consider sequences from this GTDB taxa (e.g.
    p\_\_Planctomycetota, or

<!-- -->

**\--appraise-sequence-identity** *APPRAISE_SEQUENCE_IDENTITY*

  Minimum sequence identity for SingleM appraise against reference
    database. e.g. 96% for Species-level or 86% Genus-level [default:
    0.96]

<!-- -->

**\--min-sequence-coverage** *MIN_SEQUENCE_COVERAGE*

  Minimum combined coverage for sequence inclusion [default: 10]

<!-- -->

**\--single-assembly**

  Skip appraise to discover samples to differential abundance binning.
    Forces \--num-coassembly-samples and \--max-coassembly-samples to 1
    and sets \--max- coassembly-size to None

<!-- -->

**\--exclude-coassemblies** *EXCLUDE_COASSEMBLIES* [*EXCLUDE_COASSEMBLIES* \...]

  List of coassemblies to exclude, space separated, in the form
    \"sample_1,sample_2\"

<!-- -->

**\--exclude-coassemblies-list** *EXCLUDE_COASSEMBLIES_LIST*

  List of coassemblies to exclude, space separated, in the form
    \"sample_1,sample_2\", newline separated

<!-- -->

**\--num-coassembly-samples** *NUM_COASSEMBLY_SAMPLES*

  Number of samples per coassembly cluster [default: 2]

<!-- -->

**\--max-coassembly-samples** *MAX_COASSEMBLY_SAMPLES*

  Upper bound for number of samples per coassembly cluster [default:
    \--num- coassembly-samples]

<!-- -->

**\--max-coassembly-size** *MAX_COASSEMBLY_SIZE*

  Maximum size (Gbp) of coassembly cluster [default: 50Gbp]

<!-- -->

**\--max-recovery-samples** *MAX_RECOVERY_SAMPLES*

  Upper bound for number of related samples to use for differential
    abundance binning [default: 20]

<!-- -->

**\--max-sample-combinations** *MAX_SAMPLE_COMBINATIONS*

  Maximum number of samples per target to consider for pooled clusters
    (reduces combinatorial explosions). Set to a smaller number if you
    have a polars idx error. [default: starts at 100, reducing by 20
    with retries]

<!-- -->

**\--abundance-weighted**

  Weight sequences by mean sample abundance when ranking clusters
    [default: False]

<!-- -->

**\--abundance-weighted-samples** *ABUNDANCE_WEIGHTED_SAMPLES* [*ABUNDANCE_WEIGHTED_SAMPLES* \...]

  Restrict sequence weighting to these samples. Remaining samples will
    still be used for coassembly [default: use all samples]

<!-- -->

**\--abundance-weighted-samples-list** *ABUNDANCE_WEIGHTED_SAMPLES_LIST*

  Restrict sequence weighting to these samples, newline separated.
    Remaining samples will still be used for coassembly [default: use
    all samples]

<!-- -->

**\--kmer-precluster** {never,large,always}

  Run kmer preclustering using unbinned window sequences as kmers.
    [default: large; perform preclustering when given \>1000 samples]

<!-- -->

**\--precluster-distances** *PRECLUSTER_DISTANCES*

  Distance file in the format of \`sourmash scripts pairwise\`. If
    provided, kmer sketching and clustering is skipped.

<!-- -->

**\--precluster-size** *PRECLUSTER_SIZE*

  \# of samples within each sample\'s precluster [default: 5 \*
    max-recovery- samples]

<!-- -->

**\--file-hierarchy** {never,large,always}

  Split sample outputs into subdirectories based on sample name.
    [default: large; use when given \>10000 samples]

<!-- -->

**\--file-hierarchy-depth** *FILE_HIERARCHY_DEPTH*

  Maximum depth of subdirectory hierarchy. [default: 3]

<!-- -->

**\--file-hierarchy-chars** *FILE_HIERARCHY_CHARS*

  Number of characters to use for subdirectory hierarchy. [default:
    4]

<!-- -->

**\--prodigal-meta**

  Use prodigal \"-p meta\" argument (for testing)

# COASSEMBLY OPTIONS

**\--assemble-unmapped**

  Only assemble reads that do not map to reference genomes

<!-- -->

**\--sra**

  Download reads from SRA (forward read argument intepreted as SRA
    IDs). Also sets \--run-qc.

<!-- -->

**\--download-limit** *DOWNLOAD_LIMIT*

  Parallel download limit [default: 3]

<!-- -->

**\--run-qc**

  Run Fastp QC on reads

<!-- -->

**\--unmapping-min-appraised** *UNMAPPING_MIN_APPRAISED*

  Minimum fraction of sequences binned to justify unmapping [default:
    0.1]

<!-- -->

**\--unmapping-max-identity** *UNMAPPING_MAX_IDENTITY*

  Maximum sequence identity of mapped sequences kept for coassembly
    [default: 99%]

<!-- -->

**\--unmapping-max-alignment** *UNMAPPING_MAX_ALIGNMENT*

  Maximum percent alignment of mapped sequences kept for coassembly
    [default: 99%]

<!-- -->

**\--run-aviary**

  Run Aviary commands for all identified coassemblies (unless specific
    coassemblies are chosen with \--coassemblies) [default: do not]

<!-- -->

**\--prior-assemblies** *PRIOR_ASSEMBLIES*

  Prior assemblies to use for Aviary recovery. tsv file with header:
    name [tab] assembly. Only possible with single-sample or update.
    [default: generate assemblies through Aviary assemble]

<!-- -->

**\--cluster-submission**

  Flag that cluster submission will occur through
    \`\--snakemake-profile\`. This sets the local threads of Aviary
    recover to 1, allowing parallel job submission [default: do not]

<!-- -->

**\--aviary-speed** {fast,comprehensive}

  Run Aviary recover in \'fast\' or \'comprehensive\' mode. Fast mode
    skips slow binners and refinement steps. [default: fast]

<!-- -->

**\--assembly-strategy** {dynamic,metaspades,megahit}

  Assembly strategy to use with Aviary. [default: dynamic; attempts
    metaspades and if fails, switches to megahit]

<!-- -->

**\--aviary-gtdbtk-db** *AVIARY_GTDBTK_DB*

  Path to GTDB-Tk database directory for Aviary. Only required if
    \--aviary-speed is set to comprehensive [default: use path from
    GTDBTK_DATA_PATH env variable]

<!-- -->

**\--aviary-checkm2-db** *AVIARY_CHECKM2_DB*

  Path to CheckM2 database directory for Aviary. [default: use path
    from CHECKM2DB env variable]

<!-- -->

**\--aviary-metabuli-db** *AVIARY_METABULI_DB*

  Path to MetaBuli database directory for Aviary, specifically for
    TaxVAMB. [default: use path from METABULI_DB_PATH env variable]

<!-- -->

**\--aviary-snakemake-profile** *AVIARY_SNAKEMAKE_PROFILE*

  Snakemake profile (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cli.html#profiles).
    Can be used to submit rules as jobs to cluster engine (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cluster.html).
    [default: same as \`\--snakemake-profile\`]

<!-- -->

**\--aviary-assemble-cores** *AVIARY_ASSEMBLE_CORES*

  Maximum number of cores for Aviary assemble to use. [default: 64]

<!-- -->

**\--aviary-assemble-memory** *AVIARY_ASSEMBLE_MEMORY*

  Maximum amount of memory for Aviary assemble to use (Gigabytes).
    [default: 500]

<!-- -->

**\--aviary-recover-cores** *AVIARY_RECOVER_CORES*

  Maximum number of cores for Aviary recover to use. [default: 32]

<!-- -->

**\--aviary-recover-memory** *AVIARY_RECOVER_MEMORY*

  Maximum amount of memory for Aviary recover to use (Gigabytes).
    [default: 250]

<!-- -->

**\--aviary-extra-binners** [{maxbin,maxbin2,concoct,comebin,taxvamb} \...]

  Optional list of extra binning algorithms to run. Can be any
    combination of: maxbin, maxbin2, concoct, comebin, taxvamb

<!-- -->

**\--aviary-skip-binners** [{rosella,semibin,metabat1,metabat2,metabat,vamb} \...]

  Optional list of binning algorithms to skip. Can be any combination
    of: rosella, semibin, metabat1, metabat2, metabat, vamb. Note that
    specifying

<!-- -->

**\--aviary-request-gpu**

  Request GPU resources for certain binners in Aviary recovery
    [default: do not].

# GENERAL OPTIONS

**\--output** *OUTPUT*

  Output directory [default: .]

<!-- -->

**\--cores** *CORES*

  Maximum number of cores to use [default: 1]

<!-- -->

**\--dryrun**

  dry run workflow

<!-- -->

**\--snakemake-profile** *SNAKEMAKE_PROFILE*

  Snakemake profile (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cli.html#profiles).
    Can be used to submit rules as jobs to cluster engine (see
    https://snakemake.readthedocs.io/en/v7.32.3/executing/cluster.html).

<!-- -->

**\--local-cores** *LOCAL_CORES*

  Maximum number of cores to use on localrules when running in cluster
    mode [default: 1]

<!-- -->

**\--cluster-retries** *CLUSTER_RETRIES*

  Number of times to retry a failed job when using cluster submission
    (see \`\--snakemake-profile\`) [default: 3].

<!-- -->

**\--snakemake-args** *SNAKEMAKE_ARGS*

  Additional commands to be supplied to snakemake in the form of a
    space- prefixed single string e.g. \" \--quiet\"

<!-- -->

**\--tmp-dir** *TMP_DIR*

  Path to temporary directory. [default: no default]

# OTHER GENERAL OPTIONS

**\--debug**

  output debug information

<!-- -->

**\--version**

  output version information and quit

<!-- -->

**\--quiet**

  only output errors

<!-- -->

**\--full-help**

  print longer help message

<!-- -->

**\--full-help-roff**

  print longer help message in ROFF (manpage) format

# EXAMPLES

rerun coassemble, adding new bins to database

  **\$ binchicken iterate \--coassemble-output coassemble_dir**

rerun coassemble, adding new bins to database, providing genomes directly

  **\$ binchicken iterate \--coassemble-output coassemble_dir
    \--new-genomes new_genome_1.fna**
