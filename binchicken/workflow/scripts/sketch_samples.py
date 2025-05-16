#!/usr/bin/env python3
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
import re
import argparse

SINGLEM_OTU_TABLE_SCHEMA = {
    "gene": str,
    "sample": str,
    "sequence": str,
    "num_hits": int,
    "coverage": float,
    "taxonomy": str,
    }

# Example shell directive for Snakemake:
# shell:
# """
# python3 binchicken/workflow/scripts/sketch_samples.py \
#   --unbinned {input.unbinned} \
#   --sketch {output.sketch} \
#   --taxa-of-interest {params.taxa_of_interest} \
#   --threads {threads} \
#   --log {log}
# """

def process_groups(groups, output_path):
    with SaveSignaturesToLocation(output_path) as save_sigs:
        for group in groups:
            sample, sequences = group
            mh = MinHash(n=0, ksize=60, scaled=1, track_abundance=False)
            for seq in sequences:
                mh.add_sequence(seq.replace("-", "A").replace("N", "A"))
            signature = SourmashSignature(mh, name=sample)
            save_sigs.add(signature)

def processing(unbinned_path, output_path, taxa_of_interest=None, threads=1, samples_per_group=1000):
    output_dir = os.path.dirname(output_path)

    with open(unbinned_path) as f:
        logging.info(f"Reading unbinned OTU table from {unbinned_path}")
        if taxa_of_interest:
            logging.info(f"Filtering for taxa of interest: {taxa_of_interest}")

        logging.info("Generating sketches in separate threads")
        with ProcessPoolExecutor(max_workers=threads) as executor:
            current_sample = ""
            current_sequences = []
            i = 0
            group_subset = []
            futures = []
            for line in f:
                line = line.strip()
                if line == "gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy":
                    continue

                sample = line.split("\t")[1]
                sequence = line.split("\t")[2]

                if taxa_of_interest:
                    taxonomy = line.split("\t")[5]
                    if not re.search(taxa_of_interest, taxonomy):
                        continue

                if not current_sample:
                    current_sample = sample

                if sample != current_sample:
                    group_subset.append((current_sample, current_sequences))
                    current_sample = sample
                    current_sequences = [sequence]

                    if len(group_subset) == samples_per_group:
                        output_subpath = os.path.join(output_dir, f"signatures_thread_{i}.sig")
                        future = executor.submit(process_groups, group_subset, output_subpath)
                        futures.append(future)
                        i += 1
                        group_subset = []
                else:
                    current_sequences.append(sequence)

            # Submit any remaining groups that are smaller than samples_per_group
            if sample == current_sample:
                group_subset.append((current_sample, current_sequences))

            if group_subset:
                output_subpath = os.path.join(output_dir, f"signatures_thread_{i}.sig")
                future = executor.submit(process_groups, group_subset, output_subpath)
                futures.append(future)

            for future in futures:
                future.result()

    logging.info("Concatenating sketches")
    extern.run(f"sourmash sig cat {os.path.join(output_dir, 'signatures_thread_*.sig')} -o {output_path}")

    logging.info("Done")
    return output_path

def main():
    parser = argparse.ArgumentParser(description="Sketch samples pipeline script.")
    parser.add_argument("--unbinned", required=True, help="Path to unbinned input file")
    parser.add_argument("--sketch", required=True, help="Path to output sketch file")
    parser.add_argument("--taxa-of-interest", nargs='?', const=None, default=None, help="Taxa of interest string")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars and sketching")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if args.log:
        logging.basicConfig(
            filename=args.log,
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y/%m/%d %I:%M:%S %p'
        )
    logging.info(f"Polars using {str(pl.thread_pool_size())} threads")

    processing(
        args.unbinned,
        args.sketch,
        taxa_of_interest=args.taxa_of_interest,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
