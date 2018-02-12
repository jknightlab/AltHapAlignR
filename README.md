# AltHapAlignR

* More accurate RNA-seq analysis that leverages knowledge of haplotype sequence and structure
* An approach for utilising knowledge of alternate reference haplotypes to generate gene and haplotype level estimates of transcript abundance. 




**Table of contents:**

* [Installation](#installation)
* [Mapping sequences](#mapping sequences to human genome references)
* [Usage](#usage)
* [Contact](#contact)



## Output of AltHapAlignR

### A graphic overview of haplotype prediction
Synthetic heterozygote data with the PGF and COX (1:1 ratio) haplotypes in the MHC region shown for this figure.

![](./img/output.png)


Haplotype prediction and the mapping rates (left panel). These are shown for each classical HLA gene (ordered on the y-axis according to their genomic position) with respect to each of the eight known haplotypes and presented as a heat map. Numbers in each cell are mapping rates. Predicted haplotypes are highlighted with a red border. Empty cells represent genes that are not annotated in the given haplotype. Combined mapping rates from the predicted haplotypes (middle panel). Each mapping rate in the first column is the read counts of the gene in the predicted haplotype(s) divided by the total read count of the gene across all haplotypes. Mismatching mapping rates of predicted haplotypes are in the second column. Pink and grey colors are genes predicted as heterozygous and homozygous respectively. Gene counts (right panel). Bar plots show the raw read counts for each gene.  



## Installation

Before installing 'AltHapAlignR', we need to set up python environment. 

### The python script has a few dependencies:

* [pybam](https://github.com/JohnLonginotto/pybam): "Pure Python" -but
  fast- library to read BAM files. 


```
pip install https://github.com/muffato/pybam/zipball/master
```
or
```
git clone https://github.com/JohnLonginotto/pybam    # Not available in PyPI
export PYTHONPATH=where_pybam_is${PYTHONPATH:+:$PYTHONPATH}

```
* [intervaltree](https://pypi.python.org/pypi/intervaltree): "Pure Python"
  library that implements [interval trees](https://en.wikipedia.org/wiki/Interval_tree)
```
pip install intervaltree
```
  
* [quicksect](https://pypi.python.org/pypi/quicksect): C/Python library
  that implements [interval trees](https://en.wikipedia.org/wiki/Interval_tree)
  too but is about 4x faster than `intervaltree`. Note that its
  installation may require [Cython](https://pypi.python.org/pypi/Cython)
  and a compiler (e.g. gcc) setup.
```
pip install cython      # if required
pip install quicksect
```
Only one of the last two is needed, `quicksect` being the preferred
option for performance reasons.

There are several ways of bringing them in, the easiest being with `pip`.

Note that you may want to first setup a [virtualenv](https://virtualenv.pypa.io)
before installing the dependencies, to ensure your environment is clean and
self-contained. For instance:

```sh
# if no virtualenv, 
pip install virtualenv

# Where the files are going to be stored
ALTHAPALIGN_VENV=$PWD/althapalign_virtualenv
# To create a "virtualenv" (only the first time)
virtualenv $ALTHAPALIGN_VENV
# To start using the "virtualenv"
source $ALTHAPALIGN_VENV/bin/activate

# install python modules

# To stop using it, once finished
deactivate
```

_pybam_ is unfortunately not available in PyPI, so you will have to
download it from their [github repository](https://github.com/JohnLonginotto/pybam).
The script assumes that the repository is checked out alongside, but
you can use any other location as long as it is configured in the script
itself, cf the line
```
import sys
sys.path.append('./pybam/')
```

### Installing 'AltHapAlignR'

```R

# packages to install for using AltHapAlignR :

install.packages("ggplot2", "data.table", "dplyr", "plyr", "gplots", "grid", "gridExtra", "igraph", "reshape2", "doParallel", "foreach" , "sqldf")

source("https://bioconductor.org/biocLite.R")
biocLite( c("Biostrings", "GenomicFeatures", "GenomicAlignments", "IRanges", "GenomicRanges", "Rsamtools", "rtracklayer", "") )


# install and load the 'devtools' package
install.packages("devtools")
library(devtools)
devtools::install_github('jknightlab/AltHapAlignR')

```



## Mapping sequences to human genome references


### 1. Building index reference sequences

Depending on mapper, indexing comments are different. This is an example of indexing sequences using bowtie2. 


####  building indexes of genome sequence (PGF haplotype)
`bowtie2-build hg38.genome.fa hg38.genome`
  
####   building indexes of the other haplotypes 
```bash
bowtie2-build hg38.mhc_apd.fa hg38.mhc_apd
bowtie2-build hg38.mhc_cox.fa hg38.mhc_cox
bowtie2-build hg38.mhc_dbb.fa hg38.mhc_dbb
bowtie2-build hg38.mhc_mann.fa hg38.mhc_mann
bowtie2-build hg38.mhc_mcf.fa hg38.mhc_mcf
bowtie2-build hg38.mhc_qbl.fa hg38.mhc_qbl
bowtie2-build hg38.mhc_ssto.fa hg38.mhc_ssto
```


### 2. Mapping reads

A bash script './inst/scripts/mapping2theMHCRef.sh' included in the package 




## Usage


```R
library("AltHapAlignR")

```


### Inputs

* GTF file that covers all the haplotypes
* BAM files (1 per region). They *must* be sorted by *query name* with
  [Picard](http://broadinstitute.github.io/picard/). *Do not use samtools*
  to sort the files as it [does not follow the lexicographic
  order](https://github.com/samtools/hts-specs/issues/5).



 
### Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/jknightlab/AltHapAlignR/issues>
* compose a friendly e-mail to: <wl@well.ox.ac.uk>



#### Pleases submit bugs and suggestions

This package is still under development. If you have features you would like to have added, please submit your suggestions and bug-reports at: <https://github.com/jknightlab/AltHapAlignR/issues>
