# Cockatoo

## Installation options

### Install from pip

Install latest release via pip.

```bash
pip install cockatoo-genome
```

### Install from source

Create conda env from `cockatoo.yml` and install from source.

```bash
git clone https://github.com/AroneyS/cockatoo.git
cd cockatoo
conda env create -f cockatoo.yml
conda activate cockatoo
pip install -e .
```

## Cockatoo coassemble

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).
Creates graph with samples as nodes and the number of overlapping sequences provided by SingleM.
The taxa of the considered sequences can be filtered to target a specific taxon (e.g. the phylum Planctomycetota).
The graph is clustered using the Girvan-Newman algorithm to provide sample groupings.
Aviary assemble/recover commands are generated based on proposed coassemblies.
Optionally, reads can be mapped to the matched bins with only unmapped reads being assembled.

```bash
# Example: cluster reads into proposed coassemblies based on unbinned sequences
cockatoo coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: cluster reads into proposed coassemblies based on unbinned sequences and coassemble only unbinned reads
cockatoo coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --assemble-unmapped

# Example: cluster reads into proposed coassemblies based on unbinned sequences from a specific taxa
cockatoo coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --taxa-of-interest "p__Planctomycetota"

# Example: find relevant samples for differential coverage binning (no coassembly)
cockatoo coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --single-assembly
```

## Cockatoo evaluate

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

```bash
# Example: evaluate a completed coassembly
cockatoo evaluate --coassemble-output coassemble_dir --aviary-outputs coassembly_0_dir ...
```

## Cockatoo unmap

Applies unmapping to a previous Cockatoo coassemble run, generating unmapped reads files and Aviary commands.

```bash
# Example: generate unmapped reads and commands for completed coassembly
cockatoo unmap --coassemble-output coassemble_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...
```
