---
name: binchicken
description: Use when choosing samples for metagenomic coassembly or co-binning, or running Bin Chicken (binchicken single/coassemble/iterate/update/evaluate/build) to recover novel or low-abundance MAGs from many metagenomes, including target-taxon recovery, SRA downloads, long reads and Aviary assembly/recovery.
---
# Bin Chicken

Bin Chicken suggests groups of samples for coassembly and co-binning so as to
recover genomes missing from a reference set. It profiles reads with SingleM,
appraises the marker gene windows against reference genomes (`--genomes` or
genomes recovered earlier), then clusters samples that share unbinned windows
with enough combined coverage (default ≥10X) for assembly. Assembly and
recovery can then run through Aviary.

Full docs: https://aroneys.github.io/binchicken (source in `docs/`).
Per-subcommand options: `binchicken <subcommand> --full-help`.

## Subcommands

| Subcommand   | Purpose |
|--------------|---------|
| `build`      | Create the subprocess pixi environments and set database env vars (run once) |
| `single`     | Single-sample assembly, choosing co-binning samples for each sample |
| `coassemble` | Suggest (and optionally run) coassemblies plus co-binning samples |
| `iterate`    | Rerun `coassemble`, adding bins from a previous run as references and excluding previous coassemblies |
| `update`     | Further processing of a previous run: `--run-aviary`, `--coassemblies`, `--assemble-unmapped`, `--sra` |
| `evaluate`   | Compare the targets recovered by the new genomes against what was predicted |

## Setup

```bash
binchicken build \
  --singlem-metapackage /metapackage/dir \
  --checkm2-db /checkm2/db/CheckM2_database \
  --gtdbtk-db /gtdb/release/dir \     # only for --aviary-speed comprehensive
  --metabuli-db /metabuli/db/dir      # only for --aviary-extra-binners taxvamb
```

The databases can also be set through `SINGLEM_METAPACKAGE_PATH`, `CHECKM2DB`,
`GTDBTK_DATA_PATH` and `METABULI_DB_PATH`. In a source checkout, symlink them
into `db/` instead (expected names are in `admin/set_env_vars.sh`).

## Typical workflow

```bash
# 1. Single-sample assembly with co-binning
binchicken single --forward-list fwd.txt --reverse-list rev.txt \
  --run-aviary --cores 64 --output bc_single

# 2. 2-sample coassemblies targeting windows missing from the step 1 bins
binchicken iterate --coassemble-output bc_single \
  --run-aviary --cores 64 --output bc_coassemble_2

# 3. 3-sample coassemblies
binchicken iterate --coassemble-output bc_coassemble_2 \
  --num-coassembly-samples 3 \
  --run-aviary --cores 64 --output bc_coassemble_3
```

You can also start straight from `coassemble` with your own reference genomes
(`--genomes`/`--genomes-list`). Run first without `--run-aviary` to look at the
suggestions, then rerun the same command with `--run-aviary` added, or use
`update --run-aviary --coassemblies coassembly_0 ...` to run a subset.

## Inputs

- Paired reads: `--forward`/`--reverse` (or `--forward-list`/`--reverse-list`).
  Sample names come from `*.1.fq`, `*_1.fq` or `*_R1.fq`. Reads are sorted by
  name before pairing unless `--no-sample-sort` is set.
- SRA: pass run accessions to `--forward` and add `--sra` (this also turns on
  fastp QC through `--run-qc`).
- Long reads: `--long-reads` (one file per sample). When the name matches a
  short-read sample, the long reads are merged with it. Otherwise the sample is
  long-read-only. Set the platform with `--long-read-type` (default `ont`).
  Other options are `--sra-long-reads` and `--short-long-read-pairs` (a TSV with
  columns `sample`, `long_reads`).
- Precomputed intermediates skip steps: `--sample-singlem[-dir]`,
  `--sample-query[-dir]`, `--sample-read-size`, `--genome-transcripts`,
  `--genome-singlem`, `--precluster-distances`.

## Key options

- Cluster size: `--num-coassembly-samples` / `--max-coassembly-samples`
  (default 2; max defaults to num). `--max-coassembly-size` (default 50 Gbp).
- Co-binning: `--max-recovery-samples` (default 20).
- Novelty: `--appraise-sequence-identity` (default 0.96, roughly species level;
  0.86 is roughly genus level). `--min-sequence-coverage` (default 10).
- Taxon targeting: `--taxa-of-interest` and `--exclude-taxa` take GTDB taxon
  strings or a regex such as `"p__Bacillota|p__Bacteroidota"`. The exclusion is
  applied after `--taxa-of-interest`.
- Restricting samples: `--coassembly-samples` (this is a sample list, not a
  count), `--anchor-samples`, `--exclude-coassemblies "s1,s2"`.
- Ranking: `--abundance-weighted` (experimental), optionally with
  `--abundance-weighted-samples`.
- Scale: `--kmer-precluster {never,large,always}` (preclusters automatically
  above 1000 samples), `--precluster-size` (default 5 × max-recovery-samples),
  `--file-hierarchy` (automatic above 10000 samples).
- Unmapping: `--assemble-unmapped` assembles only the reads that do not map to
  the matched references. This lowers the cost but tends to recover fewer novel
  genomes.
- `iterate`/`evaluate` bin filtering: `--min-completeness 70`,
  `--max-contamination 10` (CheckM2). `--new-genomes` sets the bins directly.

## Aviary / cluster execution

- `--run-aviary` runs assembly and recovery. Assembly uses
  `--assembly-strategy dynamic` by default: metaSPAdes at 32 CPU/250 GB, then
  64 CPU/500 GB, then MEGAHIT. Recovery runs with `--aviary-speed fast` by
  default. `comprehensive` needs GTDB-Tk.
- On HPC, add `--snakemake-profile <profile> --cluster-submission --local-cores 64 --cores 64`.
  Assemble jobs go to the cluster. Recover runs locally and Aviary handles its
  own submission (`--aviary-snakemake-profile`, which defaults to the same
  profile).
- Resources: `--aviary-assemble-cores/-memory` (64/500),
  `--aviary-recover-cores/-memory` (32/250). Failed jobs are retried with
  `--retries` (default 3), each retry with more walltime or resources.
- Without `--run-aviary`, Aviary commands are written to `coassemble/commands/`.
- Use `--dryrun` to check what will run before a long run.

## Outputs (under `--output`)

- `coassemble/target/elusive_clusters.tsv`: the suggested clusters, with
  columns `samples`, `length`, `total_targets`, `total_size`,
  `recover_samples` and `coassembly` (e.g. `coassembly_0`).
- `coassemble/target/targets.tsv`: the target windows, for checking why a
  sample was or was not chosen.
- `coassemble/summary.tsv`: the cluster summary.
- `coassemble/appraise/{binned,unbinned}.otu_table.tsv`: the SingleM appraise
  results that `iterate` reuses.
- `coassemble/coassemble/<coassembly>/{assemble,recover}`: the Aviary outputs.
- `recovered/`: the combined bins (`bins/*.fna`) with `quality_report.tsv`,
  `bin_info.tsv` and `bin_provenance.tsv`.
- `previous_bins/` (`iterate` only): symlinks to the previous iteration's
  passing bins, used as references when they are taken from Aviary outputs
  rather than from the prior run's `recovered/`.
- `evaluate/evaluate/summary_stats.tsv`: the `evaluate` summary.

## Caveats

- Novelty is defined against the SingleM metapackage and the reference genomes
  you supply. Record the metapackage and GTDB release used, because
  `--taxa-of-interest` matches against the metapackage's taxonomy.
- The 10X coverage threshold is a prediction made from marker windows. It does
  not guarantee the genome will be recovered, so use `evaluate` to check.
- Matching uses identical windows so that near-relatives do not form chimeric
  bins. Even so, larger coassemblies raise the risk of strain heterogeneity, so
  check contamination before using the bins.
- `iterate` keeps only bins that pass the completeness/contamination
  thresholds. Bins that fail do not count as references, so their lineages are
  targeted again.
- Polars error `maximum length reached`: lower `--max-sample-combinations`
  (default 100), but keep it ≥ `--max-recovery-samples`.
- Out of memory in `target_elusive`/`cluster_graph`: use
  `--kmer-precluster always`, a smaller `--precluster-size`, or narrow the
  taxa.
- `Too many open files`: run `ulimit -n 10000` before running.

## Development (source checkout)

- `pixi run postinstall` once (with `detached-environments = false` in
  `binchicken/.pixi/config.toml`), then `pixi run binchicken ...`.
- Tests: `pixi run run-a-test <name>` for a single test, then
  `pixi run run-quick-tests`, then `pixi run run-tests-at-cmr` for the expensive
  run. Tests go in the matching `test/test_<subcommand>.py`, and expensive ones
  in `test/test_manual.py`.
- Document new options in `docs/` (the `docs/tools/*.md` files are generated
  from the preludes plus `--full-help`).
