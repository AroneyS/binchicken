#################################
### collect_reference_bins.py ###
#################################
# Author: Samuel Aroney

import pandas as pd
import numpy as np
import extern

def trimmed_mean(data, trim=0.1):
    cut = int(np.floor(len(data) * trim))
    if cut == 0:
        return np.mean(data)
    else:
        a = sorted(data)
        return np.mean(a[cut:-cut])

appraise_binned = pd.read_csv(snakemake.input.appraise_binned, sep="\t")
appraise_binned["sample"] = appraise_binned["sample"].str.replace(".1$", "", regex=True)
appraise_binned = appraise_binned[appraise_binned["sample"] == snakemake.params.sample]
appraise_binned["found_in"] = appraise_binned["found_in"].str.split(",")
appraise_binned = appraise_binned.explode("found_in")

trimmed_binned = (appraise_binned.groupby(["gene", "found_in"])["coverage"]
    .sum()
    .reset_index()
    .pivot(index="gene", columns="found_in", values="coverage")
    .reset_index()
    .melt(id_vars="gene")
    .fillna(0)
    .groupby("found_in")["value"]
    .apply(trimmed_mean)
    .reset_index()
    )

reference_bins = set(trimmed_binned[trimmed_binned["value"] > 0]["found_in"].to_list())

if len(reference_bins) == 0:
    print(f"Warning: No reference bins found for {snakemake.wildcards.read}")
    cmd = f"touch {snakemake.output}"
    extern.run(cmd)
else:
    for bin in reference_bins:
        cmd = f"cat {snakemake.params.genomes[bin]} >> {snakemake.output}"
        extern.run(cmd)
