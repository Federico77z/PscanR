pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms,type = 4)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = FALSE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  
  do.call(PSMatrixList, pfms)
}

ps_get_results_table <- function(pfms)
{
  
}