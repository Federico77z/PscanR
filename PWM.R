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

J2016_pwm <- toPWM(J2016, type = "prob", pseudocounts = 0.1)
J2020_pwm <- toPWM(J2020, type = "prob", pseudocounts = 0.1)

pwm <- Matrix(J2020_pwm$MA0818.1)

hits <- matchPWM(pwm, prom_seq_unique$NR_046018.2, min.score = 0, with.score = TRUE)
