# nascentrna

[![Build Status](https://api.travis-ci.org/coregenomics/nascentrna.svg)](https://travis-ci.org/coregenomics/nascentrna)
[![Test Coverage](https://codecov.io/gh/coregenomics/nascentrna/branch/master/graphs/badge.svg)](https://codecov.io/gh/coregenomics/nascentrna)
[![License: GPL-3+](https://img.shields.io/badge/license-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

The goal of `nascentrna` is to make it easy to narrow in on interesting regions
of genome activity measured by nascent RNA sequencing
by classifying these active regions in relation to known, annotated genes
using a decision tree.


## Installation

You can install the latest version of nascentrna
from [GitHub](https://github.com/coregenomics/nascentrna) with:

``` r
# install.packages("devtools")
devtools::install_github("coregenomics/nascentrna")
```
## Example

Below is a basic example of running the `annotate()` function
to find interesting gene activity.

``` r
library(nascentrna)   # annotate()
library(rtracklayer)  # import.bed()
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Obtain transcription start sites derived from experimental nascent RNA data.
file_tss <- system.file(package = "dREG.HD", "extdata", "k562.chr21.predictions.bed")
if (file_tss == "") {
  file_tss <- "https://raw.githubusercontent.com/coregenomics/dREG.HD/master/dREG.HD/inst/extdata/k562.chr21.predictions.bed"
}
tss <- import.bed(file_tss)

# Fetch annotated genes near our TSS data.
genes_all <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_long <- subsetByOverlaps(genes_all, range(tss * 0.9))
genes <- genes_long[width(genes) < 1e6]

# Find interesting activity by proximity to genes.
annotate(tss, genes)
```
