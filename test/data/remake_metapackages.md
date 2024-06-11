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
