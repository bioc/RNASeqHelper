# RNASeqHelper for R / Bioconductor

A package to perform a full RNA-seq analysis quickly and easily.

Perform a full DESeq2 analysis for your RNA-seq data, generating
colourful Volcano and Kmeans-clustered Heatmaps, along with
time-series gene plots for genes of interest. QC-metrics that such
as PCA validation are built in, and the heatmaps generated show
z-score scaled expression between contrasts for contrasted samples
and all. Time series plots show gene trends on normalised,
corrected, and scaled data, for varying cluster correlation
thresholds.

## Installation

Install the package from Bioconductor or Gitlab, ensuring correct
dependencies.

#### From Bioconductor

```r
BiocManager::install("RNASeqHelper")
```

#### From Bioconda

```bash
micromamba install -c bioconda -c conda-forge bioconductor-rnaseqhelper
```

#### From Gitlab

```r
library(remotes)
remotes::install_github("mtekman/RNASeqHelper",
                        repos = BiocManager::repositories())
```


## Getting Started

Load the library

```r
library(RNASeqHelper)
```

Follow the vignette to learn more about how to get started with this package.
