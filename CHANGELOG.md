# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- `iterate` now falls back to `--aviary-outputs`/`--new-genomes` when the prior run's `recovered/` dir has no bins passing the quality thresholds, instead of failing with "Programming error: no bins to evaluate"
- `--assemble-unmapped` with `--file-hierarchy` (default for >=10,000 samples) collected no reference genomes, so no reads were unmapped; without `--run-qc` it failed with a `KeyError`
- `sketch_samples.py` sketched the header of appraise tables (with `found_in`) as a sample named `sample`, dropped the last kept sample when `--taxa-of-interest`/`--exclude-taxa` filtered the final rows, and merged stale per-chunk signature files from earlier runs
- `evaluate --cluster` matched coassemblies by substring (`coassembly_1` also took bins from `coassembly_10`+), and failed for single-assembly bin names
- `evaluate` now applies `--min-completeness`/`--max-contamination` to bins taken from a prior run's `recovered/` dir (previously all bins were evaluated)
- Read size counting (`read_size.csv`) now fails on corrupt/truncated read files instead of silently undercounting
- `query_processing.py` overwrites rather than appends to existing outputs
- `build --download-databases` failed on a fresh install, as the downloads ran in pixi environments (`pixi run --as-is`) before they were installed; environments are now created first
- SRA downloads delivered as a single (interleaved) file were not found when the output path contained `_1` (e.g. some temporary directories), failing with "Missing downloads"
- `evaluate` crashed (`ColumnNotFoundError: "match"`) when evaluating bins from a prior run's `recovered/` dir, as their names (`binchicken_<coassembly>.<n>`) were parsed by splitting on `-`
- `evaluate` assigned bins to the wrong coassembly when coassembly/sample names contain `-` (e.g. Aviary bins of `sample-a` and `sample-b` were merged) or start with `it<N>_`; bin names are now resolved against the known coassembly names
- `sketch_samples.py` now fails if a sample's rows are not contiguous in the unbinned table, instead of silently sketching it twice from partial sequences
- Coassembly cluster choice and recovery sample choice are now deterministic on ties (broken by `samples_hash` and sample name respectively)
- Docs: `iterate` examples used `--coassembly-samples 3` (a sample list) instead of `--num-coassembly-samples 3`; demo pointed to `recovered_bins/` instead of `recovered/bins/`

### Changed

- Smaller cluster resource requests: memory is no longer 8 GB per thread for every rule. Rules whose memory scales with input (`sketch_samples`, `distance_samples`, `abundance_weighting`, `target_weighting`, `update_appraise`) request 10x input size within 16 GB to their previous request, and others use fixed requests (e.g. `map_reads` 64 GB, `qc_reads` 16 GB, `count_bp_*`/`filter_bam_files`/`bam_to_fastq` 8 GB with 8 threads). `target_elusive` and `cluster_graph` are unchanged. Memory still doubles on each retry
- List-writing rules run locally instead of being submitted as cluster jobs, and the per-sample SingleM OTU table rename is now part of `singlem pipe` jobs (removing one cluster job per sample)
- `coverm filter` now uses its `--threads` for output compression
- Genomes are piped through SingleM in batches of 100 (one cluster job per batch rather than one per genome; new genomes in update/iterate runs are batched separately). Per-genome tables are only moved into place once the whole batch has finished. Rows within each genome's OTU table are now sorted, as SingleM's within-gene row order varied between runs
- `collect_genomes` scans the appraise tables lazily, reading only the sample's rows into memory (~4x lower peak memory per sample job)
- Faster preclustering in `target_elusive`: the top precluster samples are chosen from the sparse distance rows directly instead of densifying each row (O(n^2) overall; ~10x faster at 20k samples, 95 s at 1M), and targets are re-scanned per chunk from a temporary Parquet copy filtered to the chunk's samples. Choice of zero-similarity filler samples for samples with few neighbours is now deterministic
- `iterate` now links previous-iteration bins taken from Aviary outputs into `previous_bins/` (was `recovered_bins/`), to avoid confusion with the `recovered/` output dir

## [0.16.0] - 2026-11-24

### Added

- `--exclude-taxa` option to exclude sequences matching GTDB taxa (regex), the inverse of `--taxa-of-interest`

### Fixed

- Shell-quote `--taxa-of-interest` in workflow rules so regexes containing `|` work

### Changed

- Subworkflow environments are now run with `pixi run --as-is` instead of `--frozen`, so runs no longer write to (or contend over) already-installed environments, e.g. read-only installs. Environments are no longer installed lazily on first use: run `binchicken build` once after installing (`--skip-aviary-envs` if Aviary is not needed). CI now does this before running tests

## [0.15.1] - 2026-11-04

### Added
- Support for coassemble/single/iterate runs consisting entirely of long-read samples (`--long-reads`/`--long-reads-list`/`--sra-long-reads`/`--sra-long-reads-list`/`--short-long-read-pairs`), without requiring any `--forward`/`--forward-list`/`--sra` short reads. Long-read file sizes now also count towards `total_size`/`--max-coassembly-size`, both for long-read-only samples and for samples with both short and long reads (summed)

### Changed
- Minimum Aviary version bumped to 0.13.4, which fixes hybrid (short + long read) assembly failing in short-read Racon polishing with "empty sequences set!" when short-read headers carry a description field (e.g. SRA-downloaded reads)
- Capped the pixi dependency to <0.82: newer pixi fails every rule with "the environment directory has moved" in git checkouts, where the repo-root .pixi/pixi.toml are symlinks (prefix-dev/pixi#7335)

## [0.15.0] - 2026-09-03

### Added
- `--long-reads`/`--long-reads-list` and `--long-read-type` arguments for supplying long-read (e.g. Nanopore) sequences alongside short reads. Long reads are piped through SingleM and merged into the matched short-read sample's marker gene profile, contributing to coassembly clustering decisions, then are combined with the matched short reads for Aviary assembly and recovery. Long-read samples without a matching short-read sample are also supported, contributing to clustering and being assembled/recovered with Aviary using long reads alone
- `--sra-long-reads`/`--sra-long-reads-list` arguments for downloading long-read-only samples directly from SRA/ENA by run accession, and `--short-long-read-pairs` (a TSV file of sample/long_reads pairs, each value a local file or SRA/ENA accession) for explicitly matching a long read to an existing `--forward`/`--reverse`/`--sra` sample, bypassing the automatic filename-based matching used by `--long-reads`

### Changed
- Minimum Aviary version bumped to 0.13.3, which updates CoverM to 0.8.0 and fixes long-read coverage calculation failing with "Cannot continue without minimap2" during Aviary recovery
- Updated SingleM to v0.21.4, fixing a DIAMOND hanging bug (#223). Thanks @megan-a-wallace for reporting

## [0.14.1] - 2026-07-30

### Fixed
- Fix binchicken module not found error

## [0.14.0] - 2026-07-28

### Added
- `quickbin` to `--aviary-extra-binners` choices
- GTDB-Tk database downloading support
- `--immediate-submit` argument for immediate cluster submission mode

### Changed
- Updated Aviary to v0.13.1, ensuring internal pixi ≥ 0.72.1
- Updated CoverM to v0.8.0
- Updated SingleM to v0.21.3, restricting diamond to ≤2.2.1 for compatibility
- Updated Snakemake to v9.22
- Updated fastp to v1.3.3
- Updated Kingfisher to v0.5.0
- Pinned sra-tools to 3.2.* for compatibility
- Upgraded pixi.lock to v7 format (requires pixi ≥ 0.70)

### Fixed
- Sort read names before matching (#199). Thanks @magicprotoss for reporting
- Fix various issues stemming from CheckM1 (#193, #208). Thanks @liaohu1231 and @michoug for reporting
- Fix Aviary bin info changes due to dropping CheckM1 (#210)
- Fix single-sample target counts being incorrect (#212)
- Use proper conda packages for dependecies to avoid git clone (#196). Thanks @Asa12138 for reporting
- Fix `is_in` Categorical/Int64 crash in `cluster_graph` and `target_elusive` when sample names are purely numeric (#221). Thanks @AmaliT for reporting

## [0.13.5] - 2025-09-26

### Added
- `check_binchicken_completion.py` script for reporting Aviary assembly/recovery completion status
- `--dry-run` as alternative syntax for `--dryrun`

### Changed
- More accessible outputs including centralised bins folder and improved logging
- **Breaking**: Renamed `--cluster-retries` to `--retries`

## [0.13.4] - 2025-08-29

### Fixed
- Fix extern error type

## [0.13.3] - 2025-08-26

### Fixed
- Fixes to binchicken build environment variable setup

## [0.13.2] - 2025-08-22

### Added
- Enhanced error handling for snakemake rules with improved logging

### Fixed
- Fix build without arguments
- Minor bug fixes

## [0.13.1] - 2025-08-21

### Changed
- Set environment variables in main pixi environment
- Update Aviary to fix prepare_binning_files_split bug

## [0.13.0] - 2025-08-14

### Added
- Split preclustering into chunks to improve memory usage

### Changed
- Switch to pixi for dependency environment handling
- Enhanced retries for cluster graph to reduce combinatorial explosions
- Improved Bin Chicken vs SingleM sample naming handling

## [0.12.11] - 2025-05-19

### Fixed
- Fix PyPI publishing

## [0.12.10] - 2025-05-19

### Fixed
- Fix PyPI publishing

## [0.12.9] - 2025-05-19

### Fixed
- Fix PyPI publishing

## [0.12.8] - 2025-05-19

### Fixed
- Fix PyPI publishing

## [0.12.7] - 2025-05-19

### Added
- Conda integration for SRA and local read handling
- Anchor samples implementation
- SRA-specific arguments

### Changed
- Performance and memory usage improvements for >750,000 samples

## [0.12.6] - 2024-12-02

### Added
- Allowed provided assemblies as input
- Regex support for `--taxa-of-interest` arguments

### Changed
- Updated default appraise settings
- Bumped Aviary dependency

### Fixed
- Enumeration and suffix removal bugs

## [0.12.5] - 2024-09-06

### Added
- Separate subcommand for single-sample assembly

## [0.12.4] - 2024-08-30

### Changed
- Release to match docker image release

## [0.12.3] - 2024-08-30

### Fixed
- Minor fixes

## [0.12.2] - 2024-08-30

### Added
- Coassemblies list argument
- Database downloading support

### Changed
- Upgraded Polars to v1.x

### Fixed
- Abundance weighting calculation
- Update function bugs

## [0.12.1] - 2024-07-25

### Fixed
- Fix pip packaging

## [0.12.0] - 2024-07-24

### Added
- K-mer preclustering enabling scaling beyond 250,000 samples
- Abundance weighting biasing coassembly selection toward unrecovered lineages

## [0.11.0] - 2024-06-10

### Changed
- Updated SingleM to v0.18
- Removed GTDBtk from Aviary fast mode

### Fixed
- Fixed suffix trimming
- Fixed nontargets evaluation

## [0.10.5] - 2024-05-10

### Changed
- Dependency updates
- Aviary updates

### Fixed
- BAM naming bugs
- tmpdir corrections
- Evaluate hardening
- Coverage calculation fixes

## [0.10.4] - 2024-04-08

### Added
- Doctave-based documentation
- Sample skip notifications

### Changed
- Improved help messages
- Enhanced snakemake retry behaviour

### Fixed
- COASSEMBLY_SAMPLES sync issues

## [0.10.3] - 2024-01-15

### Added
- Links and citation

## [0.10.2] - 2024-01-15

### Fixed
- Fix PyPI upload

## [0.10.1] - 2024-01-11

### Added
- New Zenodo release

### Changed
- Increased memory for rules with retries

## [0.10.0] - 2023-12-21

### Added
- Aviary build integration
- 50 Gbp maximum assembly size

### Changed
- Project renamed to binchicken
- Restructured arguments for database paths
- Removed `--no-genomes` flag (no-genomes mode is now default)

## [0.9.6] - 2023-11-20

### Added
- Coassembly sample restriction input
- SingleM metapackage environment variable

### Changed
- Aviary environment updates

### Fixed
- Fix single assembly maximum recovery samples

## [0.9.5] - 2023-11-07

### Added
- Suffix options

### Changed
- Updated SingleM
- Improved QC handling
- Updated Polars to v0.19

### Fixed
- Fix SRA split downloads
- Prevent Aviary recover output deletion

## [0.9.4] - 2023-09-28

### Added
- Separated read downloading step

### Changed
- Switch to cores argument for Snakemake
- Switch to Fastp for QC

### Removed
- Removed bbduk QC filtering

## [0.9.3] - 2023-09-20

### Added
- Coassemble argument
- Filtered samples

### Changed
- Updated kingfisher to v0.3.0

### Fixed
- Fix contig name deduplication

## [0.9.2] - 2023-09-08

### Changed
- Enhanced Snakemake cluster support

## [0.9.1] - 2023-08-30

### Fixed
- Fix `_1` being used instead of `.1` in read file naming

## [0.9.0] - 2023-08-29

### Added
- Build subcommand

### Changed
- Switch to greedy algorithm
- Upgraded to Aviary v0.7.0

### Fixed
- Conda environment fixes
- Memory issue fixes

## [0.8.0] - 2023-05-31

### Added
- SRA interleaved read handling

### Changed
- Renamed `update` to `unmap`
- Switch to aligned fraction unmapping (99% default threshold)

## [0.7.2] - 2023-05-18

### Added
- taxa-of-interest filtering in query processing

## [0.7.1] - 2023-05-17

### Fixed
- Fix PyPI upload

## [0.7.0] - 2023-05-17

### Changed
- Updated unmapping maximum identity default to 99%
- Accounting for nontarget unbinned sequences

### Removed
- Removed S3.18.EIF_2_alpha sequences

## [0.6.0] - 2023-04-06

### Added
- Coassembly selection argument
- No-genomes option for coassembly from scratch

### Changed
- Migrated from pandas to polars

## [0.5.1] - 2023-02-23

### Added
- SRA-specific commands with download and local Aviary execution options

## [0.5.0] - 2023-02-10

### Added
- Initial release (formerly Ibis)
