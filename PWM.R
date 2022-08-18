library("TFBSTools")
library("JASPAR2020")
library("JASPAR2018")
library("JASPAR2016")

library("JASPAR2022")

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2020 <- getMatrixSet(JASPAR2020, opts)
J2016 <- getMatrixSet(JASPAR2016, opts)
