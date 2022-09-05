# SingleM coassembly

## Coassembly description

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).
Creates graph with samples as nodes and the number of overlapping sequences provided by SingleM.
The taxa of the considered sequences can be filtered to target a specific taxon (e.g. the phylum Planctomycetota).
The graph is clustered using the Girvan-Newman algorithm to provide sample groupings.

* Input: Sample reads, reference genomes
* Parameters: target taxa, coassembly sample limit, SingleM metapackage
* Output: suggested coassembly combinations

### Coassembly steps

* Create conda environment: `conda env create -f env/coassembly.yml`.
* Activate environment: `conda activate coassembly`.
* Create config with inputs and parameters e.g. `config/coassembly.yaml`.
* Run pipeline: `snakemake -s coassembly.smk --configfile config/coassembly.yaml --cores 64 --use-conda`.

## Evaluate description

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

* Input: Coassembly output, Aviary output
* Parameters: CheckM version, SingleM metapackage
* Output: suggested coassembly combinations

### Evaluate steps

* Create conda environment: `conda env create -f env/coassembly.yml`.
* Activate environment: `conda activate coassembly`.
* Create config with inputs and parameters e.g. `config/evaluate.yaml`.
* Run pipeline: `snakemake -s evaluate.smk --configfile config/evaluate.yaml --cores 64 --use-conda`.
