ps_build_bg <- function(x, pfms)
{
  checks(x, pfms)
  
  x <- BiocGenerics::unique(x)
  
  #vapply(pfms, FUN = TFBSTools::name, FUN.VALUE = character(1))
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  pfms <- lapply(pfms, FUN = ps_scan, x)
}