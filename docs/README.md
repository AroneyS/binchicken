Bin chicken
=============

Bin chicken is a tool that performs targeted recovery of low abundance metagenome assembled genomes through intelligent coassembly.

![Bin chicken workflow](/workflow.png)

## Example workflow

```bash
# Assemble and recover from each sample individually
# 20 samples used for differential abundance binning
binchicken coassemble \
  --forward-list samples_forward.txt --reverse-list samples_reverse.txt \
  --run-aviary --single-assembly \
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
