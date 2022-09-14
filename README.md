# Cockatoo

## Install from source

```bash
git clone https://github.com/AroneyS/cockatoo.git
cd cockatoo
conda env create -f cockatoo.yml
conda activate cockatoo
pip install -e .
```

## Cockatoo cluster

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).
Creates graph with samples as nodes and the number of overlapping sequences provided by SingleM.
The taxa of the considered sequences can be filtered to target a specific taxon (e.g. the phylum Planctomycetota).
The graph is clustered using the Girvan-Newman algorithm to provide sample groupings.

* Input: Sample reads, reference genomes
* Parameters: target taxa, coassembly sample limit, SingleM metapackage
* Output: suggested coassembly combinations

## Cockatoo coassemble

Snakemake pipeline to generate Aviary assemble/recover commands based on suggested coassemblies from Cluster.
Optionally, reads can be mapped to the matched bins with only unmapped reads being assembled.

* Input: Sample reads, Cluster output
* Parameters: coassemble unmapped reads flag
* Output: Aviary assemble/recover commands

## Cockatoo evaluate

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

* Input: Coassembly output, Aviary output
* Parameters: CheckM version, SingleM metapackage
* Output: suggested coassembly combinations
