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

def pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=0.1):
    appraise_binned["sample"] = appraise_binned["sample"].str.replace("\.1$", "", regex=True)
    appraise_binned = appraise_binned[appraise_binned["sample"] == sample]

    appraise_unbinned["sample"] = appraise_unbinned["sample"].str.replace("\.1$", "", regex=True)
    appraise_unbinned = appraise_unbinned[appraise_unbinned["sample"] == sample]

    num_binned = sum(appraise_binned["num_hits"].to_list())
    num_unbinned = sum(appraise_unbinned["num_hits"].to_list())
    perc_binned = num_binned / (num_binned + num_unbinned)

    if perc_binned < MIN_APPRAISED:
        return set()

    appraise_binned["found_in"] = appraise_binned["found_in"].str.split(",")
    appraise_binned = appraise_binned.explode("found_in")
    appraise_binned["found_in"] = appraise_binned["found_in"].str.replace("_protein$", "", regex=True)

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
    return reference_bins

if __name__ == "__main__":
    binned_path = snakemake.input.appraise_binned
    unbinned_path = snakemake.input.appraise_unbinned
    genomes = snakemake.params.genomes
    sample = snakemake.params.sample
    MIN_APPRAISED = snakemake.params.min_appraised
    sample_read = snakemake.wildcards.read
    output_path = snakemake.output

    appraise_binned = pd.read_csv(binned_path, sep="\t")
    appraise_unbinned = pd.read_csv(unbinned_path, sep="\t")

    reference_bins = pipeline(appraise_binned, appraise_unbinned, sample, MIN_APPRAISED=MIN_APPRAISED)

    if len(reference_bins) == 0:
        print(f"Warning: No reference bins found for {sample_read}")
        cmd = f"touch {output_path}"
        extern.run(cmd)
    else:
        for bin in reference_bins:
            cmd = f"cat {genomes[bin]} >> {output_path}"
            extern.run(cmd)