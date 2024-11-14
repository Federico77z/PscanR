#' .ps_checks
#' 
#' Validates input for motif and background sequence analysis.
#'
#' @param x Can be a `DNAStringSet` object, a character string that indicates the file path, or 
#' a `data.frame`.
#' @param pfms A list of position frequency matrices representing transcription factors motif,
#' obtained from the JASPAR database. It is expected to be a `PSMatrixList` object or a `PFMatriList` one.
#' @param type. Integer value from 1 to 4, indicating the type of check to be performed. 
#' If `type == 4`, `pfms` is expected to be a `PSMatrixList` object.
#' If `type == 3`, `x` is expected to be a `data.frame`.
#' If `type == 2`, `x` is expected to be a character string indicating a file path. 
#' If `type == 1` or `type == 4`, `x` is expected to be a `DNAStringSet` object. 
#' 
#' @details
#' The function calls `.ps_required_packages()` to ensure that the required packages are loaded. 
#' Specific checks are performed based on `type` value: 
#' * `Type == 4`: verifies that `pfms` is a `PSMatrixList` object. Computes the average 
#' and standard deviation for the background of each motif in pfms. 
#' If any values are missing (NA), it gives a warning message with the name and ID of the affected motif.
#' * `Type == 2`: if `x` is a file path (character string), checks if the file is accessible
#' and search for the `[SHORT TFBS MATRIX]` header in the first line.
#' * `Type == 1` or `Type == 4`: checks if `x` is a `DNAStringSet` object.
#' * `Type == 3`: checks if `x` is a `data.frame` containing the BG_SIZE, BG_MEANS, and BG_STDEV
#' columns. It also checks if all the values in the columns are numeric.
#' 
#' @return None. The function is used for input validation and will stop with a warning/error 
#' message if the checks are not met.
#'
#' @examples
.ps_checks <- function(x, pfms, type)
{
 # .ps_required_packages()
  
  if(type == 4 && !is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")
  
  if(type == 4)
  {
    nabg <- vapply(pfms, ps_bg_avg, FUN.VALUE = numeric(length = 1L))
    nabg2 <- vapply(pfms, ps_bg_std_dev, FUN.VALUE = numeric(length = 1L))
    
    nabg <- is.na(nabg)
    nabg2 <- is.na(nabg2)
    
    nabg <- nabg | nabg2
    
    if(any(nabg))
      warning(paste("\nNo Pscan background for", name(pfms)[nabg], ID(pfms)[nabg]))
  }
  
  if(!is(pfms, "PFMatrixList") && !is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PFMatrixList or PSMatrixList class")  
  
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
  else if(!is(x, "DNAStringSet") && (type == 1 || type == 4))
    stop("x is not an object of DNAStringSet class")
  
  else if(is.data.frame(x) && type == 3) {
    req_cols <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
    
    if(!all(req_cols %in% colnames(x)))
      stop("x does not contain required columns \"BG_SIZE\", \"BG_MEAN\", \"BG_STDEV\"")
    
    if(!all(is.numeric(x$BG_SIZE), is.numeric(x$BG_MEAN), is.numeric(x$BG_STDEV)))
      stop("Required columns of x must be of numeric type")
  }
}
#' .ps_checks2()
#' 
#' Provides additional checks on `pfms` parameter and file accessibility for writing. 
#'
#' @param pfms `PSMatrixList` object. If it is not, an error is rised.
#' @param ... Character string representing a file path. When provided, the function checks if it writable.
#' 
#' @return None, it gives errors if the checks are not met.
#'
#' @examples
.ps_checks2 <- function(pfms, ...)
{
  .ps_required_packages()

  if(!is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")  
  
 if(is.character(file))
   if(file.access(x, mode = 2) != 0)
     stop(paste("Cannot write to file path:", x))
   
}

.ps_required_packages <- function()
{
  require("Biostrings")
  require("TFBSTools")
  require("BiocParallel")
  require("BSDA")
}
