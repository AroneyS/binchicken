---
title: installation
---

Installation
========

There are several ways to install Bin Chicken

## Install from Bioconda

Install latest release via bioconda.

```bash
conda create -n binchicken -c bioconda -c conda-forge binchicken
```

## Install via Singularity / Apptainer

Install container from dockerhub.

```bash
singularity pull docker:://samuelaroney/binchicken:0.12.1
singularity run binchicken_0.12.1.sif binchicken -h
```

If your data and desired output are in the current working directory,
Bin Chicken `coassemble` can be run like so:

```bash
singularity run -B $(pwd) binchicken_0.10.5.sif binchicken coassemble \
    --forward $(pwd)/reads_1.1.fq ... \
    --reverse $(pwd)/reads_1.2.fq ...
```

Note: Bin Chicken `build` is unnecessary for this method since the conda
environments, the SingleM metapackage and CheckM2 database are included
in the container.

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
