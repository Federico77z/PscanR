ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  #pfms <- lapply(pfms, FUN = ps_scan, x)
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = TRUE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
}

ps_build_bg_from_file <- function(file, pfms)
{
  .ps_checks(file, pfms)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  short.matrix <- read.table(file, header = FALSE, row.names = 1, skip = 1)
  
  colnames(short.matrix)[1:3] <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")

  if(length(pfms) != nrow(short.matrix))
    warning("Mismatch between number of PFMs in PFMatrixList object and file table")
  
  pfms <- lapply(pfms, FUN = .ps_bg_from_file, short.matrix)
}