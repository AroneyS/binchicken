#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate aviary-v0.4.3

BASENAME=$(echo $0 | sed "s=.*/==" | sed "s/\.sh//")
CPUS=32
MEMORY=250
NCHUNKS=1
HOURS=168

while getopts 'i:o:a:c:m:n:s:t:' flag; do
  case "${flag}" in
    i) INPUT_COMMAND_LIST="${OPTARG}" ;;
    o) OUTPUT_DIR="${OPTARG}" ;;
    a) ASSEMBLE_OUTPUT="${OPTARG}" ;;
    c) CPUS="${OPTARG}" ;;
    m) MEMORY="${OPTARG}" ;;
    n) NCHUNKS="${OPTARG}" ;;
    s) SUFFIX="${OPTARG}" ;;
    t) HOURS="${OPTARG}" ;;
    *) echo "Usage: -i input_command_list -o output_dir (relative to wd) [-a assemble_output (absolute) -c cpus -m memory -n num_chunks -s suffix -t time_hours]"
       exit 1 ;;
  esac
done

BASENAME=${BASENAME}$SUFFIX
mkdir -p $OUTPUT_DIR/logs
COMMAND_FILE=$OUTPUT_DIR/logs/${BASENAME}_complete

export OUTPUT_DIR
export ASSEMBLE_OUTPUT
export CPUS
export MEMORY

envsubst < $INPUT_COMMAND_LIST > $COMMAND_FILE

mqsub --name $BASENAME --command-file $COMMAND_FILE --chunk-num $NCHUNKS --mem $MEMORY --cpus $CPUS --hours $HOURS --run-tmp-dir \
  &> $OUTPUT_DIR/logs/${BASENAME}.log
