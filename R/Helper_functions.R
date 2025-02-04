#' @keywords internal
#' @importFrom TFBSTools PFMatrixList
#' 
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
      warning(paste("\nNo Pscan background for", name(pfms)[nabg], 
                    ID(pfms)[nabg]))
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
      stop("x does not contain required columns",
           "\"BG_SIZE\", \"BG_MEAN\", \"BG_STDEV\"")
    
    if(!all(is.numeric(x$BG_SIZE), is.numeric(x$BG_MEAN), 
            is.numeric(x$BG_STDEV)))
      stop("Required columns of x must be of numeric type")
  }
}

#' @keywords internal
#' @importFrom utils globalVariables
.ps_checks2 <- function(pfms, ...)
{
  #.ps_required_packages()

  if(!is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")  
  
  if(is.character(file))
    if(file.access(x, mode = 2) != 0)
      stop(paste("Cannot write to file path:", x))
   
}

.clean_sequence <- function(x){
  
  seq_widths <- Biostrings::width(x)
  ref_width <- max(Biostrings::width(x))
  
  diff_length_seq <- x[seq_widths != ref_width]
  
  if(length(diff_length_seq) != 0){
    warning(paste(length(diff_length_seq), 'sequences found with length different from the reference. Removing the following sequences:', 
                  paste(names(diff_length_seq), collapse = ", ")))
  }
  x <- x[seq_widths == ref_width]
  
  freqs <- Biostrings::alphabetFrequency(x, as.prob = FALSE)
  n_proportions <- freqs[, "N"] / ref_width
  rem_names <- names(x[n_proportions > 0.5])
  
  if(length(rem_names)> 0){
    warning(paste('Found', length(rem_names), 'sequences with more than 50% of N. Removing the following sequences:', 
                  paste(rem_names, collapse = ", ")))
  }
  x <- x[n_proportions <= 0.5]
  return(x)
}

#.ps_required_packages <- function()
#{
#  require("Biostrings")
#  require("TFBSTools")
#  require("BiocParallel")
#  require("BSDA")
#}
