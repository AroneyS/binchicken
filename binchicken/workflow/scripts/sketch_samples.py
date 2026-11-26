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
import shlex
import tempfile
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
#   --taxa-of-interest {params.taxa_of_interest:q} \
#   --exclude-taxa {params.exclude_taxa:q} \
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

def processing(unbinned_path, output_path, taxa_of_interest=None, exclude_taxa=None, threads=1, samples_per_group=1000):
    output_dir = os.path.dirname(output_path) or "."

    # Per-chunk signatures go in a fresh directory so files left by an earlier run are never merged in
    with open(unbinned_path) as f, tempfile.TemporaryDirectory(prefix="signatures_", dir=output_dir) as chunk_dir:
        logging.info(f"Reading unbinned OTU table from {unbinned_path}")
        if taxa_of_interest:
            logging.info(f"Filtering for taxa of interest: {taxa_of_interest}")
        if exclude_taxa:
            logging.info(f"Excluding taxa: {exclude_taxa}")

        # Locate columns by name: appraise tables may carry extra columns (e.g. found_in)
        header = f.readline().rstrip("\n").split("\t")
        sample_idx = header.index("sample")
        sequence_idx = header.index("sequence")
        taxonomy_idx = header.index("taxonomy")

        logging.info("Generating sketches in separate threads")
        with ProcessPoolExecutor(max_workers=threads) as executor:
            current_sample = ""
            current_sequences = []
            i = 0
            group_subset = []
            futures = []

            def submit(groups, i):
                output_subpath = os.path.join(chunk_dir, f"signatures_thread_{i}.sig")
                futures.append(executor.submit(process_groups, groups, output_subpath))

            for line in f:
                fields = line.rstrip("\n").split("\t")
                sample = fields[sample_idx]
                sequence = fields[sequence_idx]

                if taxa_of_interest or exclude_taxa:
                    taxonomy = fields[taxonomy_idx]
                    if taxa_of_interest and not re.search(taxa_of_interest, taxonomy):
                        continue
                    if exclude_taxa and re.search(exclude_taxa, taxonomy):
                        continue

                if not current_sample:
                    current_sample = sample

                if sample != current_sample:
                    group_subset.append((current_sample, current_sequences))
                    current_sample = sample
                    current_sequences = [sequence]

                    if len(group_subset) == samples_per_group:
                        submit(group_subset, i)
                        i += 1
                        group_subset = []
                else:
                    current_sequences.append(sequence)

            # Flush the last kept sample, then any remaining groups smaller than samples_per_group
            if current_sample:
                group_subset.append((current_sample, current_sequences))

            if group_subset:
                submit(group_subset, i)

            for future in futures:
                future.result()

        logging.info("Concatenating sketches")
        extern.run(f"sourmash sig cat {shlex.quote(chunk_dir)}/signatures_thread_*.sig -o {shlex.quote(output_path)}")

    logging.info("Done")
    return output_path

def main():
    parser = argparse.ArgumentParser(description="Sketch samples pipeline script.")
    parser.add_argument("--unbinned", required=True, help="Path to unbinned input file")
    parser.add_argument("--sketch", required=True, help="Path to output sketch file")
    parser.add_argument("--taxa-of-interest", nargs='?', const=None, default=None, help="Taxa of interest string")
    parser.add_argument("--exclude-taxa", nargs='?', const=None, default=None, help="Taxa to exclude string")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for Polars and sketching")
    parser.add_argument("--log", default=None, help="Log file path")
    args = parser.parse_args()

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
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
        exclude_taxa=args.exclude_taxa,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
