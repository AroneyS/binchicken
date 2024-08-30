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
