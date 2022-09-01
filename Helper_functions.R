.ps_checks <- function(x, pfms, type)
{
  require("Biostrings")
  require("TFBSTools")
  require("BiocParallel")
  
  if(!is(pfms, "PFMatrixList") && !is(pfms, "PSMatrixList"))
  {
    stop("pfms is not an object of PFMatrixList or PSMatrixList class")  
  }
  
  if(is.character(x) && type == 2)
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
  else if(!is(x, "DNAStringSet") && type == 1)
    stop("x is not an object of DNAStringSet class")
  
  else if(is.data.frame(x) && type == 3) {
    req_cols <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
    
    if(!all(req_cols %in% colnames(x)))
      stop("x does not contain required columns \"BG_SIZE\", \"BG_MEAN\", \"BG_STDEV\"")
    
    if(!all(is.numeric(x$BG_SIZE), is.numeric(x$BG_MEAN), is.numeric(x$BG_STDEV)))
      stop("Required columns of x must be of numeric type")
    }
}

.ps_checks2 <- function(x, pfms)
{
  require("Biostrings")
  require("TFBSTools")
  require("BiocParallel")
  
 # if(!is(pfms, "PSMatrix"))
#    stop("pfms is not an object of PFMatrixList class")  
  
  if(!all(vapply(pfms, is, logical(length = 1L), "PSMatrix")))
    stop("pfms is not a list of PSMatrix objects") 
  
  if(!is(x, "DNAStringSet") && type == 1)
    stop("x is not an object of DNAStringSet class")
  
}

