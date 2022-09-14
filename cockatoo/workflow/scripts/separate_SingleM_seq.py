###############################
### separate_SingleM_seq.py ###
###############################
# Author: Samuel Aroney

import pandas as pd
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

otu_table = pd.read_csv(snakemake.input.otu_table, sep="\t")

otu_table["id"] = otu_table.index
otu_table["id"] = otu_table["id"].astype(str)

for gene in otu_table["gene"].unique():
    gene_table = otu_table[otu_table["gene"] == gene]
    records = [SeqRecord(Seq(re.sub("-", "", gene_table.iloc[i]["sequence"])), id = gene_table.iloc[i]["id"]) for i in range(len(gene_table))]
    SeqIO.write(records, snakemake.params.output_dir + "/" + gene + ".fasta", 'fasta')

otu_table.to_csv(snakemake.output.id_otu_table, sep="\t", index=False)
