
Applies further processing to a previous Bin Chicken coassemble run.
Note that all coassemblies can be run by rerunning the `coassemble` command unchanged except for adding `--run-aviary`.

Any combinations of the following:

- Generating unmapped reads files (`--assemble-unmapped`)
- Running assembly/recovery for all/specific coassemblies through Aviary (`--run-aviary`, `--coassemblies`)
- Downloading SRA reads (`--sra`)

```bash
# Example: update previous run to run specific coassemblies
binchicken update --coassemble-output coassemble_dir --run-aviary \
    --coassemblies coassembly_0 ...

# Example: update previous run to perform unmapping
binchicken update --coassemble-output coassemble_dir --assemble-unmapped

# Example: update previous run to download SRA reads
# Note: requires sample names to be SRA IDs (e.g. SRA123456)
binchicken update --coassemble-output coassemble_dir --sra

# Example: update previous run to download SRA reads, perform unmapping and run specific coassemblies
binchicken update --coassemble-output coassemble_dir --sra \
    --assemble-unmapped \
    --run-aviary --coassemblies coassembly_0 ...
```
