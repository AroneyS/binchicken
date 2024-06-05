
```bash
singularity build singularity/binchicken.sif singularity/binchicken.def
singularity run singularity/binchicken.sif binchicken -h
singularity shell -C singularity/binchicken.sif

singularity run -B /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble,$(pwd) singularity/binchicken5.sif \
    binchicken coassemble \
    --assemble-unmapped \
    --forward \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/ERR599149_1.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/ERR599166_1.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR12561417_1.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR4028167_1.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR4028175_1.fastq.gz \
    --reverse \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/ERR599149_2.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/ERR599166_2.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR12561417_2.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR4028167_2.fastq.gz \
    /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/sra/SRR4028175_2.fastq.gz \
    --run-aviary --aviary-gtdbtk-db . \
    --cores 32 \
    --genomes /work/microbiome/ibis/SRA/results/benchmarking/20240129/binchicken_co195/single_sample/coassemble/coassemble/coassembly_0/recover/bins/final_bins/metabat_bins_sens.tsv.022.fna \
    --output singularity/example/singularity_local
```

Singularity container runs fine locally
May need to modify snakemake_mqsub with an additional option (`--singularity singularity/binchicken.sif`) that sends the command to execute with the prefix of singularity exec/run
Use mqsub argument `--script-shell` with `singularity run singularity/binchicken.sif`

Alternatively, use containerize to put all the conda env in one container
But this container can't include binchicken itself
