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
import extern

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

def process_groups(groups, output_path):
    with SaveSignaturesToLocation(output_path) as save_sigs:
        for group in groups:
            sample, sequences = group
            mh = MinHash(n=0, ksize=60, scaled=1, track_abundance=False)
            for seq in sequences.iter_rows():
                mh.add_sequence(seq[1].replace("-", "A").replace("N", "A"))
            signature = SourmashSignature(mh, name=sample[0])
            save_sigs.add(signature)

def processing(unbinned, output_path, threads=1):
    output_dir = os.path.dirname(output_path)

    logging.info("Grouping samples")
    groups = [g for g in unbinned.set_sorted("sample").select("sample", "sequence").group_by(["sample"])]
    threads = min(threads, len(groups))

    # Distribute groups among threads more evenly
    grouped = [[] for _ in range(threads)]
    for i, group in enumerate(groups):
        grouped[i % threads].append(group)

    del groups

    logging.info("Generating sketches in separate threads")
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i, group_subset in enumerate(grouped):
            output_subpath = os.path.join(output_dir, f"signatures_thread_{i}.sig")
            future = executor.submit(process_groups, group_subset, output_subpath)
            futures.append(future)

        for future in futures:
            future.result()

    logging.info("Concatenating sketches")
    extern.run(f"sourmash sig cat {os.path.join(output_dir, 'signatures_thread_*.sig')} -o {output_path}")

    logging.info("Done")
    return output_path

if __name__ == "__main__":
    os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
    import polars as pl

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

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

    signatures = processing(unbinned, output_path, threads=threads)
