
Evaluates the recovery of target genes by coassemblies suggested by above, finding the number of target genes present in the newly recovered genomes.
Compares the recovery by phyla and by single-copy marker gene.

```bash
# Example: evaluate a completed coassembly
binchicken evaluate --coassemble-output coassemble_dir \
    --aviary-outputs coassembly_0_dir ...

# Example: evaluate a completed coassembly by providing genomes directly
binchicken evaluate --coassemble-output coassemble_dir \
    --new-genomes genome_1.fna ... --coassembly-run coassembly_0
```

Defaults to using genomes (from the provided coassemble outputs) with at least 70% complete and at most 10% contamination as estimated by CheckM2.
