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
