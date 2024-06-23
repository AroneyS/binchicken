#########################
### sketch_samples.py ###
#########################
# Author: Samuel Aroney

import polars as pl
import os
import logging
from sourmash import MinHash, SourmashSignature
from sourmash.sourmash_args import SaveSignaturesToLocation
from concurrent.futures import ProcessPoolExecutor

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

def process_group(group):
    sample, sequences = group
    mh = MinHash(n=0, ksize=60, scaled=1, track_abundance=False)
    for seq in sequences.iter_rows():
        mh.add_sequence(seq[1].replace("-", "A").replace("N", "A"))
    return SourmashSignature(mh, name=sample[0])

def processing(unbinned, threads=1):
    logging.info("Generating sketches")
    groups = [g for g in unbinned.select("sample", "sequence").group_by(["sample"])]

    signatures = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for signature in executor.map(process_group, groups):
            signatures.append(signature)

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
    threads = snakemake.threads

    unbinned = pl.read_csv(unbinned_path, separator="\t", schema_overrides=SINGLEM_OTU_TABLE_SCHEMA)

    if TAXA_OF_INTEREST:
        logging.info(f"Filtering for taxa of interest: {TAXA_OF_INTEREST}")
        unbinned = unbinned.filter(
            pl.col("taxonomy").str.contains(TAXA_OF_INTEREST, literal=True)
        )

    signatures = processing(unbinned, threads=threads)
    logging.info("Saving sketches to file")

    with SaveSignaturesToLocation(output_path) as save_sigs:
        for sig_obj in signatures:
            save_sigs.add(sig_obj)
