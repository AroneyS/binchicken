---
title: setup
---

Environment setup
========

Bin chicken uses separate conda environments for each subprocess.
Run `binchicken build` to create those subprocess conda environments and setup environment variables.

Conda prefix is the directory you want to contain the subprocess conda environments.
SingleM metapackage is the metapackage downloaded by SingleM using `singlem data` (see <https://github.com/wwood/singlem>).

The latter databases are required only if you want to run Aviary directly using the `--run-aviary` argument.
GTDB-Tk database is the directory containing the GTDB-Tk release (see <https://github.com/Ecogenomics/GTDBTk>).
CheckM2 database is the directory containing the CheckM2 database (see <https://github.com/chklovski/CheckM2>).
These can also be downloaded automatically by Aviary using `aviary configure --download gtdb singlem checkm2` (see <https://github.com/rhysnewell/aviary>).

```bash
binchicken build \
  --conda-prefix /path/to/conda/envs/dir \
  --singlem-metapackage /metapackage/dir \
  --gtdbtk-db /gtdb/release/dir \
  --checkm2-db /checkm2/db/dir
```

Alternatively, set directory to contain subprocess conda environments and environment variables manually.
Subprocess conda environments will be created when required.

```bash
conda env config vars set SNAKEMAKE_CONDA_PREFIX="/path/to/conda/envs"
conda env config vars set CONDA_ENV_PATH="/path/to/conda/envs"
conda env config vars set SINGLEM_METAPACKAGE_PATH="/metapackage/dir"
conda env config vars set GTDBTK_DATA_PATH="/gtdb/release/dir"
conda env config vars set CHECKM2DB="/checkm2/db/dir"
```
