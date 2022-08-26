.ps_checks <- function(x, pfms)
{
  require("Biostrings")
  require("TFBSTools")
  require("BiocParallel")
  
   if(!is(x, "DNAStringSet"))
    stop("x is not an object of DNAStringSet class")
  
  if(!is(pfms, "PFMatrixList"))
  {
    stop("pfms is not an object of PFMatrixList class")  
  }
}

