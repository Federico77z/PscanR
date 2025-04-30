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
      warning(sprintf("\nNo Pscan background for %s %s", name(pfms)[nabg], 
                     ID(pfms)[nabg]))
  }
  
  if(!is(pfms, "PFMatrixList") && !is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PFMatrixList or PSMatrixList class")  
  
  if(is.character(x) && type == 2)
  {
    if(file.access(x, mode = 4) != 0)
      stop(sprintf("Cannot access file path: %s",x))
    else
    {
      first_line <- readLines(x, n = 1)
        if(first_line != "[SHORT TFBS MATRIX]")
          stop(sprintf(x, "%s does not look like a Pscan .short_matrix file"))
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
.ps_checks2 <- function(pfms, file = NULL, ...)
{
  #.ps_required_packages()

  if(!is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")  
  
  if(!is.null(file) && !is.character(file)){
    stop('file must be a character string indicating the file path')
    
    file_dir <- dirname(file)
    
    if (!file.exists(file)) {
      if (file.access(file_dir, mode = 2) != 0) {
        stop(sprintf("Can't write to directory: %s", file_dir))
        }
      } else {
        if (file.access(file, mode = 2) != 0)
          stop(sprintf("Can't write to existing file: %s", file))
      }
  }
}

#' @keywords internal
.clean_sequence <- function(x){
  
  seq_widths <- Biostrings::width(x)
  ref_width <- max(Biostrings::width(x))
  
  diff_length_seq <- x[seq_widths != ref_width]
  
  if(length(diff_length_seq) != 0){
    warning(paste(length(diff_length_seq), 'sequences found with length 
                  different from the reference. Removing the following 
                  sequences:', 
                  paste(names(diff_length_seq), collapse = ", ")))
  }
  x <- x[seq_widths == ref_width]
  
  freqs <- Biostrings::alphabetFrequency(x, as.prob = FALSE)
  n_proportions <- freqs[, "N"] / ref_width
  rem_names <- names(x[n_proportions > 0.5])
  
  if(length(rem_names)> 0){
    warning(paste('Found', length(rem_names), 'sequences with more than 50% of N. 
                  Removing the following sequences:', 
                  paste(rem_names, collapse = ", ")))
  }
  x <- x[n_proportions <= 0.5]
  return(x)
}

#' @keywords internal
.mapping_unique_names <- function(x,pfms){

  original_names <- names(x)
  
  all_sequences_ID <- setNames(character(length(x)), original_names)
  
  # duplicate elimination
  unique_x <- BiocGenerics::unique(x)
  unique_names <- names(unique_x)  

  x_char <- as.character(x)
  unique_x_char <- as.character(unique_x)

  all_sequences_ID <- setNames(unique_names[match(x_char, unique_x_char)], 
                               original_names)

  x <- .clean_sequence(x)
  removed_sequences <- setdiff(original_names, names(x))
  
  # NA corresponds to sequences eliminated by .clean_sequence()
  all_sequences_ID[removed_sequences] <- NA

  pfms@transcriptIDLegend <- all_sequences_ID
  
  return(pfms)
}

#' @keywords internal
#' @importFrom utils download.file
.download_background <- function(file, destfile = NULL) {
  
  if (is.null(destfile)) {
    destfile <- file.path(tempdir(), file)
  }
  
  URL <- paste0(
    "https://raw.githubusercontent.com/dianabetelli/PscanR_backgrounds/refs/heads/main/BG_files/",
    file
  )
  
  utils::download.file(URL, destfile, mode = "wb")
  
  return(destfile)
}

.check_seq_duplicated <- function(x){
  for(i in seq_along(x)){
    if(names(x[i]) != x[i])
      warning(sprintf('%s will be evaluated instead of %s since they have the same 
                      promoter region', x[i], names(x[i])))
  }
}

.download_available_backgrounds <- function(destfile=NULL){
  if (is.null(destfile)) {
    destfile <- file.path(tempdir(), "AvailableBG.txt")
  }
  
  URL <- 'https://raw.githubusercontent.com/dianabetelli/PscanR_backgrounds/refs/heads/main/AvailableBG.txt'
  
  utils::download.file(URL, destfile, mode = 'wb')
  
  return(destfile)
}

#.ps_required_packages <- function()
#{
#  require("Biostrings")
#  require("TFBSTools")
#  require("BiocParallel")
#  require("BSDA")
#}
