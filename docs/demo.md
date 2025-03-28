---
title: demo
---

Demo
========

This demo will guide you through the process of running Bin Chicken on a small set of samples.
The only prerequisite is that Bin Chicken is fully installed and setup as per the [installation](installation.md) and [setup](setup.md) instructions.
With 64 CPUs and at least 100GB of RAM, the demo should take around 24 hours of running time to complete.

## Single-sample assembly with co-binning

Download 4 SRA samples and run single-sample assembly with co-binning.

```bash
mkdir ~/binchicken_demo
cd ~/binchicken_demo

binchicken single \
  --forward SRR14271365 SRR14271372 ERR2281802 ERR2281803 --sra \
  --cores 16 --output binchicken_single_assembly
```

The suggested assemblies with their respective binning samples can be found at
`binchicken_single_assembly/coassemble/target/elusive_clusters.tsv`.

The actual assembly and binning can be run by adding `--run-aviary`.
This will produce metagenome-assembled genomes (MAGs) from each sample.

```bash
binchicken single \
  --forward SRR14271365 SRR14271372 ERR2281802 ERR2281803 --sra \
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

As you might expect from checking the sample sources (e.g. [SRR14271365](https://sandpiper.qut.edu.au/run/SRR14271365)),
the SRR14271365/SRR14271372 and ERR2281802/ERR2281803 pairs are suggested for coassembly.

And, since we only have 2 samples from each group with no overlap, there are no suggested 3-sample coassemblies.

## Taxonomy targeting

Bin Chicken can also be used to target recovery of genomes from a particular taxonomic group of interest.
The taxa of interest can be specified using the `--taxa-of-interest` option, and must match GTDB taxonomy (e.g. p__Planctomycetota, or 'p__Bacillota|p__Bacteroidota').
Bin Chicken will filter the marker gene database to only include sequences from the specified taxa, thus only using marker genes from these taxa to choose coassemblies.

```bash
# Choose 2-sample coassemblies containing Mucilaginibacter
binchicken iterate \
  --coassemble-output binchicken_single_assembly \
  --cores 64 --output binchicken_2_targeting \
  --taxa-of-interest "g__Mucilaginibacter"
```

Both SRR14271365 and SRR14271372 samples contain matching Mucilaginibacter OTUs, so Bin Chicken suggests coassembly with these samples.
Since this would be the same coassembly as above, we will not run assembly/recovery here (i.e. `--run-aviary` is not used).

On the other hand, the ERR2281802 and ERR2281803 samples do not contain Mucilaginibacter OTUs with sufficient coverage, Bin Chicken does not suggest any coassemblies with these samples.
You can peruse the targets at `binchicken_2_targeting/coassemble/target/targets.tsv` to confirm this manually.

## Coassembly with co-binning

Alternatively, you can start with coassembly, optionally providing your own genomes.
Provided genomes will be used to determine the novelty of marker gene sequences in the provided metagenomes.

```bash
binchicken coassemble \
  --forward SRR14271365 SRR14271372 ERR2281802 ERR2281803 --sra \
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
