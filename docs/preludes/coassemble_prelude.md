
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
  --snakemake-profile qsub --cluster-submission --local-cores 64 --cores 64
```

Important options:

- Minimum and maximum cluster sizes can be specified (`--num-coassembly-samples` and `--max-coassembly-samples`, both default to 2)
- Maximum number of recovery samples for differential-abundance binning can be specified (`--max-recovery-samples`, default 20)
- Genomes can be provided and matching marker genes will be excluded (`--genomes`)
- Reads can be mapped to the matched bins with only unmapped reads being assembled (`--assemble-unmapped`).
- Assembly and recovery running options:
  - Run directly through Aviary (`--run-aviary`)
  - Run Aviary commands manually (see `coassemble/commands` in output)
  - Run coassemblies with differential-abundance-binning samples with the tool of your choice (see `coassemble/target/elusive_clusters.tsv` in output)
- The taxa of the considered sequences can be filtered to target a specific taxon (e.g. `--taxa-of-interest "p__Planctomycetota"`).
- Differential-abundance binning samples for single-assembly can also be found (`--single-assembly`)

Paired end reads of form \*.1.fq, \*_1.fq and \*_R1.fq, where \* represents the sample name are automatically detected and matched to their basename.
Most intermediate files can be provided to skip intermediate steps (e.g. SingleM otu tables, read sizes or genome transcripts; see `binchicken coassemble --full-help`).

## Abundance weighting (experimental)

By default, coassemblies are ranked by the number of feasibly-recovered target sequences they contain.
Instead, `--abundance-weighted` can be used to weight target sequences by their average abundance across samples.
This prioritises recovery of the most abundant lineages.
The samples for which abundances are calculated can be restricted using `--abundance-weighted-samples`.

## Kmer preclustering

Clustering groups of more than 1000 samples quickly leads to memory issues due to combinatorics.
Kmer preclustering can be used (default if >1000 samples are provided, or use `--kmer-precluster always`) to reduce the number of combinations that are considered.
This greatly reduces memory usage and allows scaling up to at least 250k samples.
Kmer preclustering can be disabled with `--kmer-precluster never`.

## Polars idx errors

If you encounter the following error: `Polars' maximum length reached. Consider installing 'polars-u64-idx'.`
You can try reducing the `--max-sample-combinations` option (default 100).
Bin Chicken ignores target sequences found in more than this number of samples to prevent combinatorial explosions.
These targets are also less useful for distinguishing between clusters.
It is not recommended to set this value lower than the number of requested recovery samples (`--max-recovery-samples`).
