---
title: installation
---

Installation
========

There are several ways to install Bin Chicken

## System requirements

The purpose of Bin Chicken is to suggest sets of samples for coassembly.
The coassembly process can use upwards of 250GB of RAM and 32 cores, so we recommend running on a HPC system if you intend to assemble and recover genomes.

Bin Chicken is supported for Linux (tested on SUSE 12.5). Specific dependencies are listed in [binchicken.yml](https://github.com/AroneyS/binchicken/blob/master/binchicken.yml).
Bin Chicken uses many independent pixi environments for different tasks. These environments are listed in [pixi.toml](https://github.com/AroneyS/binchicken/blob/master/binchicken/pixi.toml).
See [setup](/setup) for environment and database setup.

## Install from Bioconda via Pixi

Create pixi.toml file:

```toml
[workspace]
channels = ["conda-forge", "bioconda"]
name = "binchicken"
platforms = ["linux-64"]

[dependencies]
binchicken = "*"
```

Create pixi environment.

```bash
pixi install

# Either run within your current environment
pixi run binchicken -h
# Or enter the environment
pixi shell
```

## Install from Bioconda via Conda

Install latest release via conda.

```bash
conda create -n binchicken -c bioconda -c conda-forge binchicken

# Activate the environment
conda activate binchicken
```

## Install from pip

Create the environment using the `binchicken.yml` file then install from pip.

```bash
conda env create -n binchicken -f binchicken.yml
conda activate binchicken
pip install binchicken
```

## Install from source

To install from source, we recommend using [pixi](https://pixi.sh/).
Create conda env from `binchicken.yml` and install from source.

```bash
git clone https://github.com/AroneyS/binchicken.git
cd binchicken
pixi run postinstall
```

Then binchicken can be run using `pixi run` (or via `pixi shell`).

```bash
pixi run binchicken --help
```

When installed this way, binchicken is installed in an "editable" way (similar to `pip install -e .`),
meaning that any changes made to binchicken source are immediately available via the `binchicken` command.
This is useful for development and debugging.

When run this way, the databases required for binchicken (e.g. `CHECKM2DB`) can be symlinked from a `db/` directory in the binchicken repository.
An activation hook then ensures that these are available when in the pixi environments.
To do this, create a `db/` directory in the binchicken repository and symlink the required databases into it.

To check the expected database symlink names, see `admin/set_env_vars.sh` in the binchicken repository.
The advantage of this approach is that locations of the databases are not tracked in the repository, since they are specific to the computing cluster of the user.
