
Snakemake pipeline to discover samples for differential abundance binning for single-sample assembly based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).

```bash
# Example: find relevant samples for differential coverage binning (no coassembly)
binchicken single --forward reads_1.1.fq ... --reverse reads_1.2.fq ...

# Example: run proposed assemblies through aviary with cluster submission
# Create snakemake profile at ~/.config/snakemake/qsub with cluster, cluster-status, cluster-cancel, etc.
# See https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
binchicken single --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --run-aviary \
  --snakemake-profile qsub --cluster-submission --local-cores 64 --cores 64
```

Important options:

- Maximum number of recovery samples for differential-abundance binning can be specified (`--max-recovery-samples`, default 20)
- Genomes can be provided and matching marker genes will be excluded (`--genomes`)
- Reads can be mapped to the matched bins with only unmapped reads being assembled (`--assemble-unmapped`).
- Assembly and recovery running options:
  - Run directly through Aviary (`--run-aviary`)
  - Run Aviary commands manually (see `coassemble/commands` in output)
  - Run assemblies with differential-abundance-binning samples with the tool of your choice (see `coassemble/target/elusive_clusters.tsv` in output)
- The taxa of the considered sequences can be filtered to target a specific taxon (e.g. `--taxa-of-interest "p__Planctomycetota"`).

Paired end reads of form \*.1.fq, \*_1.fq and \*_R1.fq, where \* represents the sample name are automatically detected and matched to their basename.
Most intermediate files can be provided to skip intermediate steps (e.g. SingleM otu tables, read sizes or genome transcripts; see `binchicken coassemble --full-help`).

## Kmer preclustering

Clustering groups of more than 1000 samples quickly leads to memory issues due to combinatorics.
Kmer preclustering can be used (default if >1000 samples are provided, or use `--kmer-precluster always`) to reduce the number of combinations that are considered.
This greatly reduces memory usage and allows scaling up to at least 250k samples.
Kmer preclustering can be disabled with `--kmer-precluster never`.
