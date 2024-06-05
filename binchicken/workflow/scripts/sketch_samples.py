#########################
### sketch_samples.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os
import logging
from sourmash import MinHash, SourmashSignature, save_signatures

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

def processing(unbinned):
    logging.info("Generating sketches")
    parent_mh = MinHash(n=0, ksize=60, scaled=1, track_abundance=False)
    hashes = []
    samples = []
    for group in unbinned.select("sample", "sequence").group_by(["sample"]):
        mh = parent_mh.copy_and_clear()
        for row in group[1].iter_rows():
            mh.add_sequence(row[1].replace("-", "A").replace("N", "A"))

        samples.append(group[0][0])
        hashes.append(mh)

    signatures = [SourmashSignature(hashes[i], name=samples[i]) for i in range(len(hashes))]
    logging.info("Done")

    return signatures

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )

    unbinned_path = snakemake.input.unbinned
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    output_path = snakemake.output.sketch

    unbinned = pl.read_csv(unbinned_path, separator="\t", schema_overrides=SINGLEM_OTU_TABLE_SCHEMA)

    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    signatures = processing(unbinned)
    logging.info("Saving sketches to file")

    with open(output_path, "w") as f:
        save_signatures(signatures, f)
