
```bash
singularity build singularity/binchicken.sif singularity/binchicken.def
singularity run singularity/binchicken.sif binchicken -h

singularity run singularity/binchicken.sif \
    binchicken coassemble \
    --forward test/data/sample_1.1.fq test/data/sample_2.1.fq test/data/sample_3.1.fq \
    --reverse test/data/sample_1.2.fq test/data/sample_2.2.fq test/data/sample_3.2.fq \
    --genomes test/data/GB_GCA_013286235.1.fna \
    --singlem-metapackage test/data/singlem_metapackage.smpkg \
    --assemble-unmapped --prodigal-meta \
    --output singularity/example
```

Singularity container runs fine locally
May need to modify snakemake_mqsub with an additional option (`--singularity singularity/binchicken.sif`) that sends the command to execute with the prefix of singularity exec/run
Use mqsub argument `--script-shell` with `singularity run singularity/binchicken.sif`

Alternatively, use containerize to put all the conda env in one container
But this container can't include binchicken itself
