Bin Chicken
=============

Bin Chicken is a tool that performs targeted recovery of low abundance metagenome assembled genomes through strategic coassembly.

It maximises recovery of novel diversity by automatically identifying sets of samples for coassembly and co-binning.
In particular, it identifies groups of samples sharing novel marker genes that are predicted to have sufficient combined coverage for assembly and recovery of their associated genomes (10X minimum).
Identical matching of sequence windows is used to reduce the risk of forming chimeric bins from near-relatives.

It is currently designed to use metagenomic data sequenced using Illumina short-read technology.

![Bin Chicken workflow](/workflow.png)

## Example workflow

```bash
# Assemble and recover from each sample individually
# 20 samples used for differential abundance binning
binchicken single \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --run-aviary \
  --cores 64 --output binchicken_single_assembly

# Assemble and recover from 2-sample coassemblies
# Prioritising samples with genomes not previously recovered
binchicken iterate \
  --coassemble-output binchicken_single_assembly \
  --run-aviary --assemble-unmapped \
  --cores 64 --output binchicken_2_coassembly

# Perform another iteration of coassembly, with 3-samples this time
binchicken iterate \
  --coassembly-samples 3 \
  --coassemble-output binchicken_2_coassembly \
  --run-aviary --assemble-unmapped \
  --cores 64 --output binchicken_3_coassembly
```

## Help

If you have any questions or need help, please [open an issue](https://github.com/AroneyS/binchicken/issues).

## License

Bin Chicken is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
The source code is available at [https://github.com/AroneyS/binchicken](https://github.com/AroneyS/binchicken).

## Citation

Samuel T. N. Aroney, Rhys J. P. Newell, Gene W. Tyson and Ben J. Woodcroft.
Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes.
bioRxiv (2024): 2024-11. https://doi.org/10.1101/2024.11.24.625082
