# The rnaseqdea R package #

## Overview ##

`rnaseqdea` provides functions for RNA-Seq differential expression analysis.
It implements wrapper functions to call differential statistical models to
perform analysis. This package also provides a unified normalisation methods
for all used modelling methods. Some plots are implemented for these
modelling.

### Wrapper functions for statistical analysis methods ###

- DESeq2
- edgeR 
- NBPSeq
- EBSeq
- NOISeq
- SAMseq
- voom+limma
- TSPM
- hypothesis testing

### Normalisation methods ###

- Normal factor: `DESeq`, `TMM`, `RLE`, `UQ`. Use for modelling. 
- Variance stabilizing transformation (VST) and Regularized log
  transformation (RLT) implemented in package `DESeq2`.  Use for the
  classification and visualisation.

### Visualisation ###
 
- MA plot
- Volcano plot
- Histogram of p-values
- PCA plot
- MDS plot
- Boxplot

## Installation from GitHub ##

This package is only in GitHub. To install, use:

```r
devtools::install_github("wanchanglin/rnaseqdea")
```

## Usage ##

See the help pages of the package for details.
