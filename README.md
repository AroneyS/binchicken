# Ibis

[<img src="ibis_logo.png" width="50%" />](ibis_logo.png)

Ibis (bin chicken) - targeted recovery of low abundance genomes through intelligent coassembly.

## Installation options

### Install from source

Create conda env from `ibis.yml` and install from source.

```bash
git clone https://github.com/AroneyS/ibis.git
cd ibis
conda env create -f ibis.yml
conda activate ibis
pip install -e .
```

Create subprocess conda environments

```bash
ibis build --conda-prefix /path/to/conda/envs
```

Alternatively, set directory to contain subprocess conda environments

```bash
conda env config vars set SNAKEMAKE_CONDA_PREFIX="/path/to/conda/envs"
```

### Install from pip

Install latest release via pip.

```bash
pip install ibis-genome
```

## Example workflow

```bash
# Assemble and recover from each sample individually with 20 samples used for differential abundance binning
ibis coassemble \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --single-assembly --no-genomes \
  --max-recovery-samples 20 \
  --cores 64 --output ibis_single_assembly

# Assemble and recover from 2-sample coassemblies, prioritising samples with genomes not previously recovered
ibis coassemble \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --genomes-list single_assembly_genomes.txt \
  --sample-singlem-dir ibis_single_assembly/coassemble/pipe --sample-read-size ibis_single_assembly/coassemble/read_size.csv
  --assemble-unmapped \
  --max-coassembly-size 50 --max-recovery-samples 20 \
  --cores 64 --output ibis_2_coassembly

# Perform another iteration of coassembly, with 3-samples this time
ibis iterate \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --genomes-list single_assembly_genomes.txt \
  --sample-singlem-dir ibis_single_assembly/coassemble/pipe --sample-read-size ibis_single_assembly/coassemble/read_size.csv
  --coassemble-output ibis_2_coassembly --aviary-outputs ibis_2_coassembly/coassemble/coassemble/coassembly_* \
  --coassembly-samples 3 \
  --assemble-unmapped \
  --max-coassembly-size 50 --max-recovery-samples 20 \
  --cores 64 --output ibis_3_coassembly
```

## Ibis coassemble

Snakemake pipeline to discover coassembly sample clusters based on co-occurrence of single-copy marker genes, excluding those genes present in reference genomes (e.g. previously recovered genomes).
The taxa of the considered sequences can be filtered to target a specific taxon (e.g. the phylum Planctomycetota).
Assembly and recovery can be run directly, or the coassemblies with differential-abudance-binning samples can be run in the tool of your choice.
Aviary assemble/recover commands are also generated based on proposed coassemblies.
Optionally, reads can be mapped to the matched bins with only unmapped reads being assembled.
Paired end reads of form reads_1.1.fq, reads_1_1.fq and reads_1_R1.fq are automatically detected and matched to their basename.

```bash
# Example: cluster reads into proposed coassemblies
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --no-genomes

# Example: cluster reads into proposed coassemblies based on unbinned sequences
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: cluster reads into proposed coassemblies based on unbinned sequences and coassemble only unbinned reads
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --assemble-unmapped

# Example: cluster reads into proposed coassemblies based on unbinned sequences from a specific taxa
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ... --taxa-of-interest "p__Planctomycetota"

# Example: find relevant samples for differential coverage binning (no coassembly)
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --single-assembly

# Example: run proposed coassemblies through aviary with cluster submission
# Create snakemake profile at ~/.config/snakemake/qsub with cluster, cluster-status, cluster-cancel, etc.
# See https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
ibis coassemble --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --no-genomes \
  --run-aviary --aviary-gtdbtk-dir /path/to/gtdbtk_db  --aviary-checkm2-dir /path/to/checkm2_db \
  --snakemake-profile qsub --cluster-retries 3 --local-cores 64 --cores 64
```

## Ibis evaluate

Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

```bash
# Example: evaluate a completed coassembly
ibis evaluate --coassemble-output coassemble_dir --aviary-outputs coassembly_0_dir ...

# Example: evaluate a completed coassembly by providing genomes directly
ibis evaluate --coassemble-output coassemble_dir --new-genomes genome_1.fna ... --coassembly-run coassembly_0
```

## Ibis iterate

Run a further iteration of coassemble, including newly recovered bins.

```bash
# Example: rerun coassemble, adding new bins to database
ibis iterate --aviary-outputs coassembly_0_dir ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: rerun coassemble, adding new bins to database, providing genomes directly
ibis iterate --new-genomes new_genome_1.fna ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: rerun coassemble, adding new bins to database, excluding previous coassembly combinations
ibis iterate --exclude-coassemblies reads_1,reads_2 --new-genomes new_genome_1.fna ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...
```

## Ibis update

Applies further processing to a previous Ibis coassemble run: downloading SRA reads, generating unmapped reads files, and/or running assembly/recovery commands.

```bash
# Example: update previous run to download SRA reads
ibis update --sra --coassemble-output coassemble_dir --forward SRA000001 ... --genomes genome_1.fna ...

# Example: update previous run to perform unmapping
ibis update --assemble-unmapped --coassemble-output coassemble_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to run specific coassemblies
ibis update --run-aviary --coassemblies coassembly_0 ... --coassemble-output coassemble_dir --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...
```
