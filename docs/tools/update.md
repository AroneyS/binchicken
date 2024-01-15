---
title: update
---

Bin chicken update
========

Applies further processing to a previous Bin chicken coassemble run.

Any combinations of the following:

- Generating unmapped reads files (`--assemble-unmapped`)
- Running assembly/recovery through Aviary (`--run-aviary`)
- Downloading SRA reads (`--sra`)

```bash
# Example: update previous run to perform unmapping
binchicken update --coassemble-output coassemble_dir --assemble-unmapped --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to run specific coassemblies
binchicken update --coassemble-output coassemble_dir --run-aviary --coassemblies coassembly_0 ... --forward reads_1.1.fq ... --reverse reads_1.2.fq ... --genomes genome_1.fna ...

# Example: update previous run to download SRA reads
binchicken update --coassemble-output coassemble_dir --sra --forward SRA000001 ... --genomes genome_1.fna ...
```
