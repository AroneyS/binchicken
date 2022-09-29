#################################
### summarise_coassemblies.py ###
#################################
# Author: Samuel Aroney

import pandas as pd

elusive_clusters = pd.read_csv(snakemake.input.elusive_clusters, sep="\t")
summary = elusive_clusters[["coassembly", "samples", "length", "total_weight", "total_targets", "total_size"]]

if snakemake.input.read_size:
    read_size = pd.read_csv(snakemake.input.read_size, names = ["sample", "read_size"])
    read_sizes = read_size.set_index("sample").to_dict()["read_size"]
    summary["unmapped_size"] = summary["samples"].apply(lambda x: sum([read_sizes[sample] for sample in x.split(",")]))

summary.to_csv(snakemake.output.summary, sep="\t", index=False)
