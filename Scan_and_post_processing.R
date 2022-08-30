pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks2(x, pfms)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = FALSE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
}