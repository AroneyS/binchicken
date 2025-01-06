---
title: demo
---

Demo
========

## Single-sample assembly with co-binning

Download 4 SRA samples and run single-sample assembly with co-binning.

```bash
binchicken single \
  --forward SRR8334323 SRR8334324 SRR6797127 SRR6797128 --sra \
  --cores 16 --output binchicken_single_assembly
```

The suggested assemblies with their respective binning samples can be found at
`binchicken_single_assembly/coassemble/target/elusive_clusters.tsv`.

The actual assembly and binning can be run by adding `--run-aviary`.
This will produce metagenome-assembled genomes (MAGs) from each sample.

```bash
binchicken single \
  --forward SRR8334323 SRR8334324 SRR6797127 SRR6797128 --sra \
  --run-aviary --cores 64 --output binchicken_single_assembly
```

## Iterative coassembly with co-binning

The outputs from the previous command can be used as inputs to perform iterative coassembly.
This means that coassembly will be prioritised for combinations of samples with the most unrecovered diversity
(i.e. the most marker genes not present in previously recovered genomes).

We can start with 2-sample coassembly, then move to 3-sample coassembly.

```bash
# Assemble and recover from 2-sample coassemblies
binchicken iterate \
  --coassemble-output binchicken_single_assembly \
  --run-aviary --cores 64 --output binchicken_2_coassembly

# Perform another iteration of coassembly, with 3-samples this time
binchicken iterate \
  --coassembly-samples 3 --coassemble-output binchicken_2_coassembly \
  --run-aviary --cores 64 --output binchicken_3_coassembly
```

As you might expect from checking the sample sources (e.g. [SRR8334323](https://sandpiper.qut.edu.au/run/SRR8334323)),
the SRR8334323/SRR8334324 and SRR6797127/SRR6797128 pairs are suggested for coassembly.

## Coassembly with co-binning

Alternatively, you can start with coassembly, optionally providing your own genomes.

```bash
binchicken coassemble \
  --forward SRR8334323 SRR8334324 SRR6797127 SRR6797128 --sra \
  --genomes-list genomes.txt \
  --run-aviary --cores 64 --output binchicken_coassembly
```

## Testing

Tests can be run if installed from source to ensure Bin Chicken is installed correctly.
The subcommand tests can take upwards of 30 minutes to complete using a single thread.
The manual tests can take multiple days to complete, using 64 threads.

```bash
# Test single subcommand
python test/test_single.py
# Test coassemble subcommand
python test/test_coassemble.py
# Test iterate subcommand
python test/test_iterate.py

# Test downloading and coassembly/recovery with Aviary. Results stored in example/test_* directories.
python test/test_manual.py
```
