.PS_PSEUDOCOUNT <- 0.01

.checks <- function(x, pwms)
{
  require("Biostrings")
  require("TFBSTools")
  
   if(!is(x, "DNAStringSet"))
    stop("x is not an object of DNAStringSet class")
  
  if(!is(pwms, "PFMatrixList"))
  {
    stop("pwms is not an object of PFMatrixList class")  
  }
}

.ps_norm_matrix <- function(pwm)
{
  
  mx <- TFBSTools::Matrix(pwm)
  
  sums <- apply(mx, 2, sum)
  
  mx <- mx / sums
  
  mx <- mx + .PS_PSEUDOCOUNT
  
  sums <- apply(mx, 2, sum)
  
  mx <- log(mx / sums)
  
  Matrix(pwm) <- mx
  
  return(pwm)
}