# Bin chicken

[<img src="binchicken_logo.png" width="50%" />](binchicken_logo.png)

Bin chicken - targeted recovery of low abundance genomes through intelligent coassembly.

## Installation options

### Install from source

Create conda env from `binchicken.yml` and install from source.

```bash
git clone https://github.com/AroneyS/binchicken.git
cd binchicken
conda env create -f binchicken.yml
conda activate binchicken
pip install -e .
```

Create subprocess conda environments and setup environment variables.

```bash
binchicken build \
  --conda-prefix /path/to/conda/envs/dir \
  --singlem-metapackage /metapackage/dir \
  --gtdbtk-db /gtdb/release/dir \
  --checkm2-db /checkm2/db/dir
```

Alternatively, set directory to contain subprocess conda environments and environment variables manually.

```bash
conda env config vars set SNAKEMAKE_CONDA_PREFIX="/path/to/conda/envs"
conda env config vars set CONDA_ENV_PATH="/path/to/conda/envs"
conda env config vars set SINGLEM_METAPACKAGE_PATH="/metapackage/dir"
conda env config vars set GTDBTK_DATA_PATH="/gtdb/release/dir"
conda env config vars set CHECKM2DB="/checkm2/db/dir"
```

### Install from pip

Install latest release via pip.

```bash
pip install binchicken
```

## Example workflow

```bash
# Assemble and recover from each sample individually with 20 samples used for differential abundance binning
binchicken coassemble \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --run-aviary --single-assembly \
  --cores 64 --output binchicken_single_assembly

# Assemble and recover from 2-sample coassemblies, prioritising samples with genomes not previously recovered
binchicken iterate \
  --coassemble-output binchicken_single_assembly \
  --run-aviary --assemble-unmapped \
  --cores 64 --output binchicken_2_coassembly

# Perform another iteration of coassembly, with 3-samples this time
binchicken iterate \
  --coassembly-samples 3 \
  --coassemble-output binchicken_2_coassembly \
  --run-aviary --assemble-unmapped \
  --cores 64 --output binchicken_3_coassembly
```

## Bin chicken coassemble

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).
The taxa of the considered sequences can be filtered to target a specific taxon (e.g. the phylum Planctomycetota).
Assembly and recovery can be run directly, or the coassemblies with differential-abudance-binning samples can be run in the tool of your choice.
Aviary assemble/recover commands are also generated based on proposed coassemblies.
Optionally, reads can be mapped to the matched bins with only unmapped reads being assembled.
Paired end reads of form reads_1.1.fq, reads_1_1.fq and reads_1_R1.fq are automatically detected and matched to their basename.

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
  --snakemake-profile qsub --cluster-retries 3 --local-cores 64 --cores 64
```

## Bin chicken evaluate

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

```bash
# Example: evaluate a completed coassembly
binchicken evaluate --coassemble-output coassemble_dir --aviary-outputs coassembly_0_dir ...

# Example: evaluate a completed coassembly by providing genomes directly
binchicken evaluate --coassemble-output coassemble_dir --new-genomes genome_1.fna ... --coassembly-run coassembly_0
```

## Bin chicken iterate

Run a further iteration of coassemble, including newly recovered bins.
Defaults to using genomes with at least 70% complete and at most 10% contamination CheckM2.
Automatically excludes previous coassemblies.

```bash
# Example: rerun coassemble, adding new bins to database
binchicken iterate --coassemble-output coassemble_dir

# Example: rerun coassemble, adding new bins to database, providing genomes directly
binchicken iterate --coassemble-output coassemble_dir --new-genomes new_genome_1.fna
```

## Bin chicken update

Applies further processing to a previous Bin chicken coassemble run: downloading SRA reads, generating unmapped reads files, and/or running assembly/recovery commands.

```bash
# Example: update previous run to download SRA reads
binchicken update --coassemble-output coassemble_dir --sra --forward SRA000001 ... --genomes genome_1.fna ...

# Example: update previous run to perform unmapping
binchicken update --coassemble-output coassemble_dir --assemble-unmapped --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to run specific coassemblies
binchicken update --coassemble-output coassemble_dir --run-aviary --coassemblies coassembly_0 ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...
```
