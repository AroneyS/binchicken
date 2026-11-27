# Remake metapackages

## base metapackge

```bash
singlem makedb \
    --otu-table test/data/singlem_otu_table.tsv \
    --db test/data/singlem.sdb

singlem metapackage \
    --singlem-packages test/data/singlem_metapackage.smpkg/S3.7.ribosomal_protein_S7.spkg \
    --no-taxon-genome-lengths \
    --nucleotide-sdb test/data/singlem.sdb \
    --metapackage test/data/singlem_metapackage2.smpkg

rm -r test/data/singlem.sdb
rm -r test/data/singlem_metapackage.smpkg
mv test/data/singlem_metapackage2.smpkg test/data/singlem_metapackage.smpkg
```

## EIF metapackage

```bash
# Create SingleM package from S3.18.EIF_2_alpha
cp ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/hmmseq/S3.18.EIF_2_alpha* test/data
# Manually edited to remove most sequences

singlem regenerate \
    --input-singlem-package ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/packages/S3.7.ribosomal_protein_S7.spkg \
    --input-sequences test/data/S3.18.EIF_2_alpha.faa \
    --input-taxonomy test/data/S3.18.EIF_2_alpha_taxonomy.tsv \
    --euk-sequences ~/m/abisko/aroneys/uniprot/uniref100_20210616/uniprot_sprot.fa \
    --euk-taxonomy ~/m/abisko/aroneys/uniprot/uniref100_20210616/uniprot_sprot_taxonomy.tsv \
    --output-singlem-package test/data/S3.18.EIF_2_alpha.spkg \
    --sequence-prefix S3.18.EIF_2_alpha~

# SingleM otu table made from otu tables at ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/transcripts/
singlem makedb \
    --otu-table test/data/singlem_otu_table_EIF.tsv \
    --db test/data/singlem_EIF.sdb

singlem metapackage \
    --singlem-packages test/data/singlem_metapackage_EIF.smpkg/S3.7.ribosomal_protein_S7.spkg test/data/singlem_metapackage_EIF.smpkg/S3.18.EIF_2_alpha.spkg \
    --no-taxon-genome-lengths \
    --nucleotide-sdb test/data/singlem_EIF.sdb \
    --metapackage test/data/singlem_metapackage_EIF2.smpkg

rm -r test/data/singlem_EIF.sdb
rm -r test/data/singlem_metapackage_EIF.smpkg
mv test/data/singlem_metapackage_EIF2.smpkg test/data/singlem_metapackage_EIF.smpkg
```

## Live SingleM tables (pipe_live)

`test/data/pipe_live/` holds real `singlem pipe` outputs for the test reads and genome, so tests not about
SingleM itself can pass them with `--sample-singlem`/`--genome-singlem` and skip piping (results are identical).
Regenerate them whenever SingleM or the test metapackage changes:

```bash
binchicken coassemble \
    --forward test/data/sample_{1,2,3}.1.fq \
    --reverse test/data/sample_{1,2,3}.2.fq \
    --genomes test/data/GB_GCA_013286235.1.fna \
    --singlem-metapackage test/data/singlem_metapackage.smpkg \
    --prodigal-meta \
    --output pipe_live_run \
    --snakemake-args "$PWD/pipe_live_run/coassemble/summarise/bins_summarised.otu_table.tsv \
        $PWD/pipe_live_run/coassemble/pipe/sample_1_read.otu_table.tsv \
        $PWD/pipe_live_run/coassemble/pipe/sample_2_read.otu_table.tsv \
        $PWD/pipe_live_run/coassemble/pipe/sample_3_read.otu_table.tsv"

cp pipe_live_run/coassemble/pipe/sample_{1,2,3}_read.otu_table.tsv test/data/pipe_live/
cp pipe_live_run/coassemble/summarise/bins_summarised.otu_table.tsv \
    test/data/pipe_live/bins_summarised_GB_GCA_013286235.1.otu_table.tsv
```
