# PONG2 – KIR Genotype Imputation & Model Training

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release](https://img.shields.io/github/v/release/NormanLabUCD/PONG2)](https://github.com/NormanLabUCD/PONG2/releases)
<!-- Add real CI badge if you set up GitHub Actions later -->
<!-- [![R-CMD-check](https://github.com/NormanLabUCD/PONG2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NormanLabUCD/PONG2/actions) -->

**PONG2** is an R package with C++ acceleration (via Rcpp) for **high-accuracy imputation** and **training** of Killer-cell Immunoglobulin-like Receptor (**KIR**) genotypes from SNP array in the KIR locus (chromosome 19).

It is optimized for population genetics, immunogenetics, and large-scale biobank studies requiring reliable KIR allele calls.

**Main CLI commands:**

- `impute` – predict KIR alleles from target PLINK files  
- `train` – build a new prediction model from reference genotypes + known KIR calls

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
  - [From GitHub (recommended for latest version)](#from-github-recommended-for-latest-version)
  - [From source tarball](#from-source-tarball)
  - [From source directory (developer mode)](#from-source-directory-developer-mode)
- [CLI Setup & Quick Start](#cli-setup--quick-start)
- [Usage](#usage)
  - [impute command](#impute-command)
  - [train command](#train-command)
- [Improving Imputation Accuracy](#improving-imputation-accuracy)
- [Input & Output Formats](#input--output-formats)
- [Dependencies & External Tools](#dependencies--external-tools)
- [Troubleshooting](#troubleshooting)
- [Development & Contributing](#development--contributing)
- [License](#license)
- [Citation](#citation)
- [Contact & Support](#contact--support)

## Overview

PONG2 enables scalable and accurate KIR genotyping by combining:

- region-specific PLINK2 preprocessing  
- optional local minimac4 imputation for missing variants  
- supervised allele prediction models tailored to the polymorphic KIR region

It supports both hg19 and hg38 assemblies and is particularly useful for studying immune response variation, HLA–KIR interactions, and disease association in diverse populations.

## Features
- Automatic handling of hg19 / hg38 coordinate differences  
- Configurable SNP missingness threshold  
- Built-in local imputation fallback (`--fill-missing`)  
- Support for external pre-imputation (e.g. Michigan Imputation Server)  
- Multi-threading via `--threads`  
- Force-run mode and missing SNP imputation strategies  
- Clean separation of preprocessing and prediction steps

## Requirements

**R version**  
≥ 4.0

**Required R packages** (loaded at runtime):
- readr
- tidyverse
- parallel

**System tools** (must be in PATH):

- PLINK 2.0 (required)
- minimac4 ≥ 4.1.6 (required only when using `--fill-missing`)
- bgzip & tabix (HTSlib – usually bundled with minimac4)

## Installation

### From GitHub (recommended for latest version)

Install the development or stable version directly from the repository:

```r
# Install remotes if needed
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install PONG2 from GitHub
# Replace 'yourusername' with your actual GitHub username
remotes::install_github("https://github.com/NormanLabUCD/PONG2")
```

### Standard install
```bash
R CMD INSTALL PONG2_1.0.0.tar.gz
```

### Custom library path
```bash
R CMD INSTALL --library="installation path" PONG2_1.0.0.tar.gz
```

### USAGE
```bash
pong2 <command> [options]
```
### Basic imputation
```bash
pong2 impute -i bfile -o output -l KIR -a hg19
```

### Imputation with local missing SNP fill-in
```bash
pong2 impute -i bfile -o output -l KIR -a hg19 --fill-missing -t 20
```

### Train a new model
```bash
pong2 train -i bfile -k kfile -o output -l KIR -a hg19 -t 20
```
### Specific detailed help
```bash
pong2 --help                # General overview + list of commands
pong2 --help impute         # Detailed help for imputation
pong2 --help train          # Detailed help for training
pong2 version               # Show version number
```

### impute command
```bash
pong2 impute [options]
```
| Flag | Description | Example Value |
| :--- | :--- | :--- |
| `-i, --input` | PLINK bed/bim/fam prefix (should cover KIR locus on chr19) | `data/my_genotypes_chr19` |
| `-o, --output` | Output directory (will be created if it doesn't exist) | `results/imputation` |
| `-l, --locus` | Target locus (currently only KIR is supported) | `KIR` |
| `-a, --assembly` | Genome build used in the data | `hg19` or `hg38` |

### Optional flags

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--filter` | `0.01` or `0.005` | KIR genotype quality filter threshold used in model training |
| `-t, --threads` | `20` | Number of CPU threads |
| `-f, --force` | `false` | Proceed even if SNP matching rate is low |
| `--fill-missing` | `false` | Impute missing SNPs locally with minimac4 |


### train command
```bash
pong2 train [options]
```
| Flag | Description |
| :--- | :--- |
| `-i, --bfile` | Reference PLINK bed/bim/fam prefix |
| `-k, --kfile` | CSV with sample IDs and KIR calls (e.g., sample KIR3DL1_h1 KIR3DL1_h2) |
| `-o, --output` | Directory to save trained model |
| `-l, --locus` | KIR |
| `-a, --assembly` | hg19 or hg38 |

### Optional flags
| Flag | Default | Description |
| :--- | :--- | :--- |
| `-t, --threads` | `20` | Number of CPU threads |
| `--filter` | `0.01` or `0.005` | KIR genotype quality filter threshold |
| `--pos` | *Optional* | Optional KIR region |

### Examples
```bash
pong2 impute -i example/chr19 -o results/run1 -l KIR3DL1 -a hg19
Pre-imputation (recommended for best accuracy)
Uses minimac4 + bundled reference panels.
pong2 impute -i example/chr19 -o results/run1 -l KIR3DL1 -a hg19 --fill-missing -t 32

pong2 impute -i example/chr19 -o results/run1 -l KIR3DL1 -a hg19 -f
pong2 train -i example/chr19 -k example/kir_calls.csv -o models/v2 -l KIR -a hg19 -t 24
```

### License

```markdown
## License

PONG2 is licensed under the **GNU General Public License v3.0** (GPL-3.0).

Full license text: [LICENSE](LICENSE) (included in the repository root)

You are free to use, modify, and distribute PONG2, provided that:

- You distribute any derivative work under the same GPL-3.0 license
- You include the original copyright notice and license text

See the [GNU GPL-3.0 page](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.

## Citation

If you use PONG2 in your research, please cite:

> Sadeeq, S. A., Leaton, L., Castelli, E., & Norman, P. (2025).  
> PONG 2.0: Allele imputation for the killer cell immunoglobulin-like receptors.  
> *Human Immunology*, 86, 111488.  
> https://doi.org/10.1016/j.humimm.2025.111488
```
### BibTeX entry:

```bibtex
@article{sadeeq2025pong,
  author  = {Sadeeq, S. A. and Leaton, L. and Castelli, E. and Norman, P.},
  title   = {PONG 2.0: Allele imputation for the killer cell immunoglobulin-like receptors},
  journal = {Human Immunology},
  volume  = {86},
  pages   = {111488},
  year    = {2025},
  doi     = {10.1016/j.humimm.2025.111488}
}
```



### Contact & Support

```markdown
## Contact & Support

- **GitHub Issues** (preferred for bug reports, feature requests, questions):  
  https://github.com/NormanLabUCD/PONG2/issues

- **Email** (for collaboration, private inquiries, or urgent support):  
  paul.norman@cuanschutz.edu

- **Lab / Institution**:  
  Norman Lab / [University of Colorado / Department of Biomedical Informatics, Immunology & Microbiology]

We aim to respond to issues within 1–3 business days.

Happy KIR analysis! 🧬
```






