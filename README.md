# AltHapAlignR
=================

**Table of contents:**

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Contact](#contact)



## Output of AltHapAlignR

A graphic overview of haplotype prediction.
Synthetic heterozygote data with the PGF and COX (1:1 ratio) haplotypes in the MHC region shown for this figure.

![](./img/output.pdf)


Haplotype prediction and the mapping rates (left panel). These are shown for each classical HLA gene (ordered on the y-axis according to their genomic position) with respect to each of the eight known haplotypes and presented as a heat map. Numbers in each cell are mapping rates. Predicted haplotypes are highlighted with a red border. Empty cells represent genes that are not annotated in the given haplotype. Combined mapping rates from the predicted haplotypes (middle panel). Each mapping rate in the first column is the read counts of the gene in the predicted haplotype(s) divided by the total read count of the gene across all haplotypes. Mismatching mapping rates of predicted haplotypes are in the second column. Pink and grey colors are genes predicted as heterozygous and homozygous respectively. Gene counts (right panel). Bar plots show the raw read counts for each gene.  


## Introduction



## Installation

```R

# packages to install for using AltHapAlignR :

install.packages("ggplot2", "data.table", "dplyr", "plyr", gplots", "grid", "gridExtra", "igraph", "reshape2", "doParallel", "foreach" , "sqldf")

source("https://bioconductor.org/biocLite.R")
biocLite( c("Biostrings", "GenomicFeatures", "GenomicAlignments", "IRanges", "GenomicRanges", "Rsamtools", "rtracklayer", "") )


# install and load the 'devtools' package
install.packages("devtools")
library(devtools)
devtools::install_github('jknightlab/AltHapAlignR')

```


## Usage

```R
library("AltHapAlignR")

```


## Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/jknightlab/AltHapAlignR/issues>
* compose a friendly e-mail to: <wl@well.ox.ac.uk>



## Pleases submit bugs and suggestions

This package is still under development. If you have features you would like to have added, please submit your suggestions and bug-reports at: <https://github.com/jknightlab/AltHapAlignR/issues>