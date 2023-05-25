###########################
### query_processing.py ###
###########################
# Author: Samuel Aroney

import os
import sys
import argparse
import logging
import gzip
import polars as pl

READS_COLUMNS={
    "read": str,
    "number": int,
    }

def is_interleaved(reads, total_count):
    """
    Return True if:
     - reads have duplicate metadata in consecutive reads
     - total_count is even
    """
    READ_METADATA_REGEX = r"[^ ]*(.*)"

    if total_count % 2 == 1:
        return False, "Odd readcount"

    processed = (
        reads
        .select(
            pl.col("number").add(1).floordiv(2).alias("index"),
            pl.col("number").mod(2).alias("rank"),
            pl.col("read").str.extract(READ_METADATA_REGEX)
            )
        .pivot(values="read", index="index", columns="rank", aggregate_function=None)
        .filter(pl.col("0") != pl.col("1"))
    )

    if processed.height > 0:
        return False, f"Consecutive reads do not match ({processed.height}/{reads.height // 2})"

    return True, "Duplicate read names in consecutive reads with even readcount"

def pipeline(fastq_reads, START_CHECK_PAIRS, END_CHECK_PAIRS):
    start_list = []
    end_list = []
    read_count = 0

    if fastq_reads.endswith(".gz"):
        f = gzip.open(fastq_reads, "rt")
    else:
        f = open(fastq_reads)

    for line_count, line in enumerate(f):
        if line_count % 4 != 0:
            continue

        read_count += 1
        if len(start_list) < START_CHECK_PAIRS * 2:
            start_list.append([line.strip(), read_count])

        end_list.append([line.strip(), read_count])
        if len(end_list) > END_CHECK_PAIRS * 2:
            end_list.pop(0)

    f.close()

    try:
        line_count
    except UnboundLocalError:
        return False, "Empty file"

    reads = pl.DataFrame(start_list + end_list, schema=READS_COLUMNS)

    return is_interleaved(reads, read_count)

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--input", help="Input fastq file (can be gzipped)")
    parser.add_argument("--output", type=argparse.FileType("w"), default=sys.stdout, help="Output file")

    parser.add_argument("--start-check-pairs", type=int, default=5, help="Number of pairs of reads to compare at the start of the file")
    parser.add_argument("--end-check-pairs", type=int, default=5, help="Number of pairs of reads to compare at the end of the file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl
    logging.info(f"Polars using {str(pl.threadpool_size())} threads")

    output_outcome, output_reason = pipeline(args.input, args.start_check_pairs, args.end_check_pairs)

    args.output.write(f"{output_outcome}\t{output_reason}\n")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
