ps_build_bg <- function(x, pwms)
{
  checks(x, pwms)
  
  x <- BiocGenerics::unique(x)
  
  #vapply(pwms, FUN = TFBSTools::name, FUN.VALUE = character(1))
  
  pwms <- lapply(pwms, FUN = .ps_norm_matrix)
  
}