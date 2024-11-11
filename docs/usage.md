---
title: usage
---

General usage
========

## General Options

- **Cores**: The `--cores` option sets the maximum number of cores to use. The default value is 1. When running in cluster mode, this restricts the number of cores per job submitted.
- **Local Cores**: The `--local-cores` option sets the maximum number of cores to use on local rules when running in cluster mode. The default value is 1.
- **Cluster Retries**: The `--cluster-retries` option sets the number of times to retry a failed job. The default value is 3. This is especially important when using `--run-aviary` to run assembly and genome recovery. Jobs are submitted with sequentially increasing walltime, and sometimes CPUs and RAM with each failure. See below.

## Memory issues

If you encounter memory issues with target_elusive or cluster_graph, you can try the following:

- **Kmer preclustering**: The `--kmer-precluster` option can be set to `always` to enable kmer preclustering. This is the default if more than 1000 samples are provided. This greatly reduces memory usage and allows scaling up to at least 250k samples on a single HPC node. Kmer preclustering can be disabled with `--kmer-precluster never`.
- **Precluster size**: The `--precluster-size` option sets the maximum number of samples to use for preclustering. The default value is 5 x the number of recovery samples. Reducing this value will reduce memory usage by reducing the number of combinations that are considered. However, this reduces the samples combinations that are considered for both coassembly and differential-abundance binning (co-binning), so may result in sub-optimal coassembly suggestions.
- **Target taxa**: The `--taxa-of-interest` option can be used to filter the taxa of the considered sequences to target a specific taxon. This can reduce memory usage.

## Aviary Assemble

- **Maximum number of cores**: The `--aviary-assemble-cores` option sets the maximum number of cores for Aviary assemble to use. The default value is 64.
- **Maximum amount of memory**: The `--aviary-assemble-memory` option sets the maximum amount of memory for Aviary assemble to use (in Gigabytes). The default value is 500.
- **Dynamic assembly strategy**: The `--assembly-strategy dynamic` (default) option means that assembly is initially attempted using metaSPAdes with 32 CPUs, 250 GB RAM. If this fails, metaSPAdes is attempted again with 64 CPUs, 500 GB RAM. If this also fails, MEGAHIT is attempted with 32 CPUs, 250 GB RAM.
- **MetaSPAdes assembly strategy**: The `--assembly-strategy metaspades` option means that assembly is attempted using metaSPAdes with 32 CPUs, 250 GB RAM, with CPUs and RAM increasing with each attempt.
- **MEGAHIT assembly strategy**: The `--assembly-strategy megahit` option means that assembly is attempted using MEGAHIT with 32 CPUs, 250 GB RAM, with CPUs and RAM increasing with each attempt.

## Aviary Recover

- **Maximum number of cores**: The `--aviary-recover-cores` option sets the maximum number of cores for Aviary recover to use. The default value is 32. If `--cluster-submission` is set, then this applies to jobs submitted by Aviary. Binners, refinement steps and characterisation are submitted as separate jobs.
- **Maximum amount of memory**: The `--aviary-recover-memory` option sets the maximum amount of memory for Aviary recover to use (in Gigabytes). The default value is 250. If `--cluster-submission` is set, then this applies to jobs submitted by Aviary.
