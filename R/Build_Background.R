#' Build background matrices in parallel
#' 
#' Generates a background probability profile matrix list using the BiocParallel
#' framework for parallel execution.
#' 
#' @param x A `DNAStringSet` object (see Biostrings package) containing the set of 
#' regulatory sequences from co-regulated or co-expressed genes (i.e. a set of 
#' gene promoters). These sequences are the target for background scanning.
#' 
#' @param pfms A list of position frequency matrices representing transcription factor motifs,
#' obtained from the JASPAR database. Should be of class `PSMatrixlist`, otherwise gets 
#' automatically converted. 
#' 
#' @param BPPARAM Parallelization parameter passed to `bplapply`. See BiocParallel.
#' 
#' @param BPOPTIONS Optional configuration settings passed to `bplapply`. See BiocParallel.
#'
#' @details A check on inputs is performed with the helper function `.ps_check`, specifiyng 
#' `type == 1`, duplicated sequences are removed from `x`, `pfms` is converted into `PSMatrixList` object,
#' and the scoring of motif in parallel is performed by `bplapply()`.
#' 
#' @return A PSMatrixList object, containing each motif matrix from `pfms`, background-scored 
#' against the sequences in `x`. 
#'
#' @export
#'
#' @examples
#' # Get promoter sequences
#' prom_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoters(txdb, upstream = 200, downstream = 50))
#'
#' # Load JASPAR motif matrices for vertebrates
#' opts <- list(collection = "CORE", tax_group = "vertebrates")
#' J2020 <- getMatrixSet(JASPAR2020, opts)
#'
#' # Generate the background-scored motif matrices
#' bg_matrices <- ps_build_bg(prom_seq, J2020, BPPARAM = MulticoreParam(24))
ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms, type = 1)
  
  x <- BiocGenerics::unique(x)
  
 # if(is(pfms, "PFMatrixList"))
    pfms <- as(pfms, "PSMatrixList")
    #pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = TRUE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  
  do.call(PSMatrixList, pfms)
}
#' Build background matrices from a file 
#' 
#' Generates a background probability profile matrix list by reading an input file.
#' This function is useful when background information are stored into an external file. 
#' 
#' @param file A character string representing the path for input file. 
#' This file contains background information for regulatory sequences. Can be generated 
#' with the `ps_build_bg` function.
#' 
#' @param pfms A list of position frequency matrices representing transcription factor motif, 
#' obtained from the JASPAR database.
#'
#' @details Validates the input with the internal function `.ps_checks`, with `type` set to 2. 
#' Reads the input file into a `data.frame` with the row names taken from the first column. 
#' Calls `ps_build_bg_from_table` with the created `data.frame` and `pfms` as input.
#' 
#' @return A `PSMatrixList` object, containing each motif matrix from `pfms`, background-scored 
#' using the values provided in `file`. See `ps_build_bg_from_table`.
#' 
#' @export
#' 
#' @examples
#' # Load a background information file
#' file_path <- "../background.txt"
#'
#' # Load JASPAR motif matrices for vertebrates
#' opts <- list(collection = "CORE", tax_group = "vertebrates")
#' J2020 <- getMatrixSet(JASPAR2020, opts)
#'
#' # Generate the background-scored motif matrices from file
#' bg_matrices <- ps_build_bg_from_file(file_path, J2020)
ps_build_bg_from_file <- function(file, pfms)
{
  .ps_checks(file, pfms, type = 2)
  
  short.matrix <- read.table(file, header = FALSE, row.names = 1, skip = 1)
  
  colnames(short.matrix)[1:3] <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
  
  pfms <- ps_build_bg_from_table(short.matrix, pfms)
}
#' Build background matrices from table
#' 
#' Generates a background probability profile matrix list from a background parameter table.
#' 
#' @param x A `data.frame` containing background parameters for each motif.
#' 
#' @param pfms A list of position frequency matrices (PFMs) representing transcription factor 
#' motifs, typically from the JASPAR database. Each element should be coercible to the 
#' `PSMatrix` class. The number of elements in `pfms` should match the number of row
#' of `x`
#' 
#' @details Calls the internal function `.ps_checks()` to validate inputs, converts each
#' elements of pfms into `PSMatrix` class, and applies the internal function `.ps_bg_from_table`
#' to each element of `pfms` with the data in `x`. 
#' A warning is issued if the number of elements in `pfms` doesn't match the number of row of `x`.
#' 
#' @return A `PSMatrixList` object containing each motif matrix from `pfms`, scored with background 
#' parameters derived from `x`.
#' 
#' @export
#'
#' @examples
#' 
ps_build_bg_from_table <- function(x, pfms)
{
  .ps_checks(x, pfms, type = 3)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  if(length(pfms) != nrow(x))
    warning("Mismatch between number of PFMs in PFMatrixList/PSMatrixList object and file table")
  
  pfms <- lapply(pfms, FUN = .ps_bg_from_table, x)
  
  do.call(PSMatrixList, pfms)
}
#' Compute background statistics for position frequency matrices
#' 
#' Generates a background statistic table (size, mean, and standard deviation)
#' from a list of position frequency matrices.
#'
#' @param pfms A list of position frequency matrices representing the transcription factor motifs.
#'
#' @details Calls the helper functions `ps_bg_size`, `ps_bg_avg`, and `ps_bg_std_dev` to compute
#' the statistics. 
#' In addiction, it also validates the input with `.ps_chechs()` function.
#'
#' @return A dataframe with one row for each PFM in `pfms`. 
#' Columns: 
#' `BG_SIZE`: an integer vector representing the background size for each PFMs. 
#' `BG_MEAN`: a numeric vector representing the mean of the background frequencies for each PFMs.
#' `BG_STDEV`: a numeric vector representing the standard deviation of the background frequencies
#' for each PFMs. 
#' 
#' @export
#'
#' @examples
ps_get_bg_table <- function(pfms)
{
  .ps_checks2(pfms)
  
  BG_SIZE <- vapply(pfms, ps_bg_size, integer(length = 1)) 
  BG_MEAN <- vapply(pfms, ps_bg_avg, numeric(length = 1))
  BG_STDEV <- vapply(pfms, ps_bg_std_dev, numeric(length = 1))
  
  data.frame(BG_SIZE, BG_MEAN, BG_STDEV, row.names = names(pfms))
}
#' Write background statistics to a file
#'
#' This function computes background statistics for each element of a list of position frequency
#' matrices, and write the result in a specified file. 
#'
#' @param pfms A list of position frequency matrices representing transcription factor motif, 
#' obtained from the JASPAR database.
#'
#' @param file A character string specifing the path where the result should be saved.
#'
#' @details It validates the input by calling the heper function `ps.checks2()`.
#' By calling ps_get_bg_table, it computes the background statistics. 
#' A header is added to the file (`[SHORT TFBS MATRIX`]). 
#'
#' @return None. It saves the given background statistics in a file.
#'
#' @export
#'
#' @examples
#' 
ps_write_bg_to_file <- function(pfms, file)
{
  .ps_checks2(pfms, file)
  
  tab <- ps_get_bg_table(pfms)
  
  write("[SHORT TFBS MATRIX]", file = file, append = FALSE)
  
  write.table(tab, file = file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
}
