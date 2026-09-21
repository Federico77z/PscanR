# JASPAR 2020 CORE vertebrate PFMs (746 matrices).
# Source and license: https://jaspar.elixir.no/ (CC BY 4.0).
# Citation: Fornes et al., https://doi.org/10.1093/nar/gkz1001.
# Rscript J2020.R OUTPUT.rds; an existing file is never overwritten.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, !file.exists(args[[1]]))
matrices <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,
    list(collection = "CORE", tax_group = "vertebrates"))
saveRDS(matrices, args[[1]], version = 3)
