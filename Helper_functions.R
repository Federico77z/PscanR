.ps_checks <- function(x, pfms)
{
  require("Biostrings")
  require("TFBSTools")
  require("BiocParallel")
  
  if(is.character(x))
  {
    if(file.access(x, mode = 4) != 0)
      stop(paste("Cannot access file path:",x))
    else
    {
      first_line <- readLines(x, n = 1)
        if(first_line != "[SHORT TFBS MATRIX]")
          stop(paste(x, " does not look like a Pscan .short_matrix file"))
    }
  }
  else if(!is(x, "DNAStringSet"))
    stop("x is not an object of DNAStringSet class")
  
  if(!is(pfms, "PFMatrixList"))
  {
    stop("pfms is not an object of PFMatrixList class")  
  }
}

