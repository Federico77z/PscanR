## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# Load required libraries 
library(PscanR)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(TFBSTools)
library(txdbmaker)
library(GenomicFeatures)

