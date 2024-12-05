#' Build background matrices in parallel
#' 
#' Generates a background probability profile matrix list using the BiocParallel
#' framework for parallel execution. It computes background scores for each 
#' motif matrix from the JASPAR database against a set of regulatory sequences 
#' (e.g., gene promoters) using the `ps_scan` function.
#' 
#' @param x A `DNAStringSet` object (see Biostrings package) containing the set 
#'   of regulatory sequences from co-regulated or co-expressed genes (i.e. a set 
#'   of gene promoters). These sequences are the target for background scanning.
#' 
#' @param pfms A list of position frequency matrices representing transcription 
#'   factor motifs, obtained from the JASPAR database. Should be of class 
#'   `PSMatrixlist`, otherwise gets automatically converted. 
#' 
#' @param BPPARAM Parallelization parameter passed to `bplapply` function from. 
#'   the `BiocParallel` package. This parameter controls how the parallel
#'   processing is executed.
#'   see `BiocParallel` package for more details.
#' 
#' @param BPOPTIONS Optional configuration settings passed to `bplapply` 
#'   function from the `BiocParallel` package. This can include additional
#'   options to control the behavior of parallel execution.
#'   See `BiocParallel` documentation for more details.
#'
#' @details 
#' A check on inputs is performed with the helper function `.ps_check`, 
#' specifiyng `type == 1`. Duplicated sequences are removed from `x` to avoid 
#' redundant computations. The motif matrices in `pfms` are converted into 
#' a `PSMatrixList` class object, 
#' The motif are background scored by the `ps_scan` function in parallel.
#' 
#' @return 
#' A PSMatrixList object, containing each motif matrix from `pfms`, 
#' background-scored against the sequences in `x`. 
#'
#' @importFrom txdbmaker makeTxDbFromUCSC 
#' @importFrom GenomicFeatures promoters
#' @importFrom Biostrings getSeq
#' @importFrom TFBSTools getMatrixSet
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import JASPAR2020
#' 
#' @examples
#' # Import GTF annotation from UCSC
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#'
#' # Generate the background-scored motif matrices
#' bg_matrices <- ps_build_bg(prom_seq, J2020, BPPARAM = BiocParallel::MulticoreParam(1))
#' bg_matrices
#' @export
ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms, type = 1)
  
  x <- BiocGenerics::unique(x)
  
 # if(is(pfms, "PFMatrixList"))
    pfms <- as(pfms, "PSMatrixList")
    #pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  pfms <- bplapply(
    pfms, 
    FUN = ps_scan, 
    x, 
    BG = TRUE, 
    BPPARAM=BPPARAM, 
    BPOPTIONS = BPOPTIONS
  )
  
  do.call(PSMatrixList, pfms)
}
#' Build background matrices from a file 
#' 
#' Generates a background probability profile matrix list by reading 
#' an input file containing background information for regulatory sequences.
#' This function is useful when background information are stored into an 
#' external file. 
#' 
#' @param file A character string representing the path for input file. 
#'   This file contains background information for regulatory sequences, 
#'   such as size, mean, and standard deviation for background distributions.
#'   The file should be in tabular format, where the first column contains 
#'   the sequence identifiers (row names) and subsequent columns contain the 
#'   background statistics (e.g., `BG_SIZE`, `BG_MEAN`, `BG_STDEV`).
#'   This file can be generated with the `ps_build_bg` function.
#' 
#' @param pfms A list of position frequency matrices representing transcription 
#'   factor motif, obtained from the JASPAR database. These matrices are 
#'   used to search motif in regulatory sequences, and will be background-scored
#'   using the information from the input file. 
#'
#' @details 
#' This function:
#' 
#' - Validates the input with the internal function `.ps_checks`, 
#'   with `type` set to 2. 
#' - Reads the input file into a `data.frame`, with row names taken 
#'   from the first column. 
#' - Calls `ps_build_bg_from_table` with the created `data.frame` 
#'   and `pfms` as input.
#' 
#' @return 
#' A `PSMatrixList` object, containing each motif matrix from `pfms`, 
#' background-scored using the values provided in `file`. 
#' See `ps_build_bg_from_table` for more details.
#' 
#' @importFrom TFBSTools getMatrixSet
#' 
#' @examples
#' # Load a background information file
#' file_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR")
#'
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#'
#' # Generate the background-scored motif matrices from file
#' bg_matrices <- ps_build_bg_from_file(file_path, J2020)
#' bg_matrices
#' @seealso \code{\link{ps_build_bg_from_table}}
#' 
#' @export
ps_build_bg_from_file <- function(file, pfms)
{
  .ps_checks(file, pfms, type = 2)
  
  short.matrix <- read.table(file, header = FALSE, row.names = 1, skip = 1)
  
  colnames(short.matrix)[seq_len(3)] <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
  
  pfms <- ps_build_bg_from_table(short.matrix, pfms)
}
#' Build background matrices from table
#' 
#' Generates a background probability profile matrix list by applying 
#' background parameters from a given table to a list of motif matrices.
#' 
#' @param x A `data.frame` containing background parameters for each motif.
#'    The number of rows in `x` should match the number of motif in `pfms`.
#' 
#' @param pfms A list of position frequency matrices (PFMs) representing 
#'   transcription factor motifs, typically from the JASPAR database. 
#'   Each element should be coercible to the `PSMatrix` class. 
#'   The number of elements in `pfms` should match the number of row of `x`
#' 
#' @details 
#' This function: 
#' 
#' - Calls the internal function `.ps_checks()` to validate inputs, with 
#'   `type == 3`.
#' - Converts each elements of pfms into `PSMatrix` class. 
#' - A warning is issued if the number of elements in `pfms` doesn't match 
#'   the number of row of `x`.
#' - Applies the internal function `.ps_bg_from_table` to each element of 
#'   `pfms` with the correspponding background data from `x`. 
#'   
#' 
#' @return 
#' A `PSMatrixList` object containing each motif matrix from `pfms`, 
#' scored with background parameters provided from `x`.
#' 
#' @importFrom TFBSTools getMatrixSet
#'
#' @examples
#' # create the `data.frame`
#' background_data <- data.frame(
#' BG_SIZE = c(500, 450, 480),
#' BG_MEAN = c(0.3, 0.25, 0.35),
#' BG_STDEV = c(0.05, 0.07, 0.06)
#' )
#' 
#' # Retrieve motif matrices for vertebrates from JASPAR2020
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' J2020_subset <- J2020[1:3] # match the number of rows in `backgound_data`
#' 
#' rownames(background_data) <- c("MA0004.1", "MA0006.1", "MA0019.1")
#' 
#' # Generate background-scored motif matrices
#' bg_matrices <- ps_build_bg_from_table(background_data, J2020_subset)
#' bg_matrices
#' @export
ps_build_bg_from_table <- function(x, pfms)
{
  .ps_checks(x, pfms, type = 3)
  
  pfms <- lapply(pfms, FUN = as, "PSMatrix")
  
  if(length(pfms) != nrow(x))
    warning(
      "Mismatch between number of PFMs in PFMatrixList/PSMatrixList object", 
      "and file table")
  
  pfms <- lapply(pfms, FUN = .ps_bg_from_table, x)
  
  do.call(PSMatrixList, pfms)
}
#' Compute background statistics for position frequency matrices
#' 
#' Generates a background statistic table (size, mean, and standard deviation)
#' from a list of position frequency matrices (PFMs) representing transcription
#' factor motifs.
#'
#' @param pfms A list of position frequency matrices representing the 
#'   transcription factor motifs. Each element should be a `PSMatrix` object 
#'   (or should be coercible to `PSMatrix`).
#'   
#' @details 
#' This function computes background statistics for each motif matrix in `pfms`.
#' Before computing the statistics, it validates the input with `.ps_chechs()` 
#' function.
#' It calls the helper functions `ps_bg_size`, `ps_bg_avg`, and `ps_bg_std_dev` 
#' to compute the statistics. 
#'
#' @return 
#' A dataframe with one row for each PFM in `pfms`.
#' 
#' Columns: 
#' 
#' - `BG_SIZE`: An integer vector representing the background size for each PFMs. 
#' - `BG_MEAN`: A numeric vector representing the mean of the background 
#'   frequencies for each PFMs.
#' - `BG_STDEV`: a numeric vector representing the standard deviation of the 
#'   background frequencies for each PFMs. 
#'
#' @importFrom TFBSTools getMatrixSet
#' 
#' @examples
#' # Retrieve motif matrices for vertebrates from JASPAR2020
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' J2020_subset <- J2020[1:3] # match the number of rows in `backgound_data`
#' 
#' # create the `data.frame`
#' background_data <- data.frame(
#' BG_SIZE = c(500, 450, 480),
#' BG_MEAN = c(0.3, 0.25, 0.35),
#' BG_STDEV = c(0.05, 0.07, 0.06)
#' )
#' 
#' rownames(background_data) <- c("MA0004.1", "MA0006.1", "MA0019.1")
#' 
#' # Generate background-scored motif matrices
#' bg_matrices <- ps_build_bg_from_table(background_data, J2020_subset)
#' bg_table <- ps_get_bg_table(bg_matrices)
#' bg_table 
#' @seealso \code{\link{ps_generics}}
#' 
#' @export
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
#' Computes background statistics for each element of a list of position 
#' frequency matrices, and write the result to a specified file. 
#'
#' @param pfms A list of position frequency matrices representing transcription 
#'   factor motif, obtained from the JASPAR database. Each element should be
#'   of `PSMatrix` class or coercible to `PSMatrix`.
#'
#' @param file A character string specifying the path to the output file where 
#'   the background statistics should be saved.
#'
#' @details 
#' It validates the input by calling the heper function `ps.checks2()`.
#' It computes the background statistics (size, mean, and standard deviation)
#' using `ps_get_bg_table()`, then writes the result to a specified file. 
#' A header is added to the file (`[SHORT TFBS MATRIX`]). 
#'
#' @return None. It saves the given background statistics to the specified file
#' in a tab-delimited format.
#' 
#' @importFrom TFBSTools getMatrixSet
#'
#' @examples
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' # File path to save the result
#' file_path <- "J2020_hg38_bg_stats.txt"
#' 
#' PSM1 <- PSMatrix(
#'   pfm = J2020[[1]],
#'   ps_bg_avg = 0.25,        
#'   ps_fg_avg = 0.5,         
#'   ps_bg_std_dev = 0.05,    
#'   ps_bg_size = 250L        
#'   )
#' PSM2 <- PSMatrix(
#'   pfm = J2020[[2]],
#'   ps_bg_avg = 0.25,        
#'   ps_fg_avg = 0.5,         
#'   ps_bg_std_dev = 0.05,    
#'   ps_bg_size = 250L        
#'   )
#'   
#' PSMatrixList_J2020 <- PSMatrixList(PSM1, PSM2)
#' 
#' ps_write_bg_to_file(PSMatrixList_J2020, file_path)
#' 
#' @seealso \code{\link{ps_get_bg_table}}
#' 
#' @export
ps_write_bg_to_file <- function(pfms, file)
{
  .ps_checks2(pfms, file)
  
  tab <- ps_get_bg_table(pfms)
  
  write("[SHORT TFBS MATRIX]", file = file, append = FALSE)
  
  write.table(
    tab, 
    file = file, 
    quote = FALSE, 
    sep = "\t", 
    row.names = TRUE,
    col.names = FALSE, 
    append = TRUE
  )
}
