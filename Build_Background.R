ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms, type = 1)
  
  x <- BiocGenerics::unique(x)
  
 # if(is(pfms, "PFMatrixList"))
    pfms <- as(pfms, "PSMatrixList")
    #pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = TRUE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  
  do.call(PSMatrixList, pfms)
}

ps_build_bg_from_file <- function(file, pfms)
{
  .ps_checks(file, pfms, type = 2)
  
  short.matrix <- read.table(file, header = FALSE, row.names = 1, skip = 1)
  
  colnames(short.matrix)[1:3] <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
  
  pfms <- ps_build_bg_from_table(short.matrix, pfms)
}

ps_build_bg_from_table <- function(x, pfms)
{
  .ps_checks(x, pfms, type = 3)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  if(length(pfms) != nrow(x))
    warning("Mismatch between number of PFMs in PFMatrixList/PSMatrixList object and file table")
  
  pfms <- lapply(pfms, FUN = .ps_bg_from_table, x)
  
  do.call(PSMatrixList, pfms)
}
