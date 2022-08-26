ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  #pfms <- lapply(pfms, FUN = ps_scan, x)
  pfms <- bplapply(pfms, FUN = ps_scan, x, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
}