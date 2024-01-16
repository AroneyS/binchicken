
Run a further iteration of coassemble, including newly recovered bins.
All [coassemble](/tools/coassemble) options are available and can be altered for the new iteration (see [example workflow](/)).

```bash
# Example: rerun coassemble, adding new bins to database
binchicken iterate --coassemble-output coassemble_dir

# Example: rerun coassemble, adding new bins to database, providing genomes directly
binchicken iterate --coassemble-output coassemble_dir --new-genomes new_genome_1.fna
```

Defaults to using genomes (from the provided coassemble outputs) with at least 70% complete and at most 10% contamination as estimated by CheckM2.
Automatically excludes previous coassemblies.
