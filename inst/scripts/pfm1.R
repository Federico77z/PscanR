# Reproduce the MA0506.1 foreground fixture from bundled human inputs.
# See fixture-provenance.md for upstream sources and historical limitations.
# Rscript pfm1.R OUTPUT.rds; no downloads and no existing file is overwritten.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, !file.exists(args[[1]]))
stopifnot(requireNamespace("PscanR", quietly = TRUE))
input <- function(name) system.file("extdata", name, package="PscanR", mustWork=TRUE)
motifs <- readRDS(input("J2020.rds"))
background <- PscanR::ps_retrieve_bg_from_file(
    input("J2020_hg38_200u_50d_UCSC.psbg.txt"), motifs)["MA0506.1"]
result <- PscanR::pscan(readRDS(input("prom_seq.rds")), background,
    BPPARAM=BiocParallel::SerialParam())
saveRDS(result[[1]], args[[1]], version=3)
