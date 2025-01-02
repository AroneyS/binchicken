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

## Install from Bioconda

Install latest release via bioconda.

```bash
conda create -n binchicken -c bioconda -c conda-forge binchicken
```

## Install from pip

Install latest release via pip.

```bash
pip install binchicken
```

## Install from source

Create conda env from `binchicken.yml` and install from source.

```bash
git clone https://github.com/AroneyS/binchicken.git
cd binchicken
conda env create -f binchicken.yml
conda activate binchicken
pip install -e .
```

### Testing and demo

Tests can be run if installed from source to ensure Bin Chicken is installed correctly.
The subcommand tests can take upwards of 30 minutes to complete using a single thread.
The manual tests can take multiple days to complete, using 64 threads.

```bash
# Test single subcommand
python test/test_single.py
# Test coassemble subcommand
python test/test_coassemble.py
# Test iterate subcommand
python test/test_iterate.py

# Test downloading and coassembly/recovery with Aviary. Results stored in example/test_* directories.
python test/test_manual.py
```
