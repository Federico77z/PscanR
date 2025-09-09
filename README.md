# Overview

**PscanR** is an R/Bioconductor package for motif enrichment analysis in regulatory DNA sequences. It is thought as a more flexible and integrable solution of *pscan* algorithm, originally developed by Zambelli et al. (2009) and aviable as webserver and C++ application (<http://www.beaconlab.it/pscan>). The algorithm performs transcription factor binding motifs enrichment analysis in promoters of co-regulated or co-expressed genes by comparing matches of a Position Weight Matrix (PWM) database in the input set against a background set composed of all the promoters from the same organism. 

PscanR is implemented as a Bioconductor-ready package leveraging core functionalities from estabilished packages such as *Biostrings*, *TFBSTools*, *BiocParallel*, and others. A key feature of our implementation is the extension of the *PFMatrixList* S4 class from TFBSTools into a new PSMatrixList class, which handles most of the data processing and manipulation tasks.

# Installation from GitHub

```{r}
devtools::install_github("Federico77z/PscanR")
```

# Vignettes

For in-depth examples and instructions, please use:

```{r}
browseVignettes('PscanR')
```
