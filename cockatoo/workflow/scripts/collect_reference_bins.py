#################################
### collect_reference_bins.py ###
#################################
# Author: Samuel Aroney

import pandas as pd
import extern

appraise_binned = pd.read_csv(snakemake.input.appraise_binned, sep="\t")
appraise_binned["found_in"] = appraise_binned["found_in"].str.split(",")
appraise_binned = appraise_binned.explode("found_in")

reference_bins = set(appraise_binned["found_in"].to_list())

if len(reference_bins) == 0:
    print(f"Warning: No reference bins found for {snakemake.wildcards.read}")
    cmd = f"touch {snakemake.output}"
    extern.run(cmd)
else:
    for bin in reference_bins:
        cmd = f"cat {snakemake.params.genomes[bin]} >> {snakemake.output}"
        extern.run(cmd)
