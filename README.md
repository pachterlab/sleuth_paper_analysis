# sleuth paper analysis

This repo contains all of the code to reproduce the results in the [sleuth preprint](http://biorxiv.org/content/early/2016/06/10/058164).


# Preliminaries

- Install [snakemake](https://bitbucket.org/johanneskoester/snakemake)
- Download and install `R` along with dependencies listed below (R dependencies section)
- Updated the `BASE` variable in `config.py` to represent the base path on your system

# Organization

The code is organized into a few different directories, each with a theme:

- `annotation`: pulls down the different annotations used and creates indeces
- `bottomly`: analysis related to the Bottomly et al. data, particular the 'self-consistency FDR' experiments
- `cuffdiff2_analysis`: analysis of the Trapnell et al. dataset to extract effect sizes from that dataset
- `geuvadis`: the bulk of the simulations, based on results from the GEUVADIS data
- `simulation_core`: dependencies for the simulations in the `geuvadis` directory
- `software`: the bulk of the software used, not including the R dependencies

# R dependencies

Install using `install.packages()`

### from CRAN

- `cowplot`
- `devtools`
- `dplyr`
- `data.table`
- `ggplot2`
- `jsonlite`
- `reshape2`
- `scales`

### from Bioconductor

First, install Bioconductor:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
```

Then, you should be able to install packages using the `biocLite()` function.

- `biomaRt`
- `BitSeq`
- `DESeq`
- `DESeq2`
- `EBSeq`
- `edgeR`
- `limma`

### from GitHub

- `[sleuth v0.28.1](https://github.com/pachterlab/sleuth/tree/bioRxiv)` fork with some modifications: `devtools::install_github('pachterlab/sleuth', ref = 'bioRxiv')`
- `mamabear v0.2`: `devtools::install_github('pimentel/mamabear', ref = 'v0.2')`

# Bug reports

Please make them in [GitHub issues](https://github.com/pachterlab/sleuth_paper_analysis/issues).
