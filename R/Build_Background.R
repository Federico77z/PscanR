#' Build background matrices in parallel
#' 
#' This function generates a background probability profile matrix list, 
#' used to assess if a transcription factor binding motif is appearing more 
#' often than by chance in the pscan() input sequences, using the BiocParallel 
#' framework for parallel execution. 
#' It computes background scores for each motif matrix, collected, for example, 
#' from the JASPAR database, by scanning them against a set of regulatory 
#' sequences (e.g., gene promoters) using the `ps_scan` function.
#' Note that the promoter region analysed may vary. 
#' 
#' @param x A `DNAStringSet` object (see Biostrings package) containing the set 
#'   of all the regulatory sequences of the organism of study 
#'   (e.g., gene promoters retrieved from a `TxDb` object). 
#'   These sequences are the target for background scanning.
#' 
#' @param pfms A `PFMatrixlist` object containing position frequency matrices 
#'   (PFM) representing transcription factor binding preferences. 
#'   Those are sourced from databases such as JASPAR.
#' 
#' @param BPPARAM Parallelization parameter passed to `bplapply` function from 
#'   the `BiocParallel` package. This parameter defines the parallel processing 
#'   settings, including the number of cores or workers that the function uses 
#'   for computation. 
#'   See `BiocParallel` package for more details.
#' 
#' @param BPOPTIONS Optional configuration settings passed to `bplapply` 
#'   function from the `BiocParallel` package. This can include additional
#'   options to control the behavior of parallel execution.
#'   See `BiocParallel` documentation for more details.
#'  
#' @param fullBG Logical. Default is FALSE. When set to TRUE, it creates a 
#'   mapping between all sequence names in the organism of study 
#'   and the corresponding names retained after applying the unique() function. 
#'   For example, if multiple identical sequences exist (e.g., ID1, ID2, ID3, 
#'   and ID4), and unique() retains only ID2, the mapping will associate each 
#'   original name with its unique counterpart (ID1 → ID2, ID2 → ID2, 
#'   ID3 → ID2, ID4 → ID2). This vector is stored in the transcriptIDLegend 
#'   slot of the PSMatrixList and helps generate a complete background 
#'   PSMatrixList, improving computational efficiency for the pscan_fullBG() 
#'   function.
#'   In addition, it retrieves for each PSM in the PSMatrixList output all the 
#'   background metrics relative to hits score, position, strand, and 
#'   oligonucleotide sequence of the regulatory sequences scanned with the PWM. 
#'     
#' @details 
#' This function validates input types and removes duplicated sequences
#' from `x` to avoid redundant computations. It also remove all the sequences 
#' with an N content above 50%.
#' The motif matrices are background scored by the `ps_scan` function in 
#' parallel.
#' The output of this function can be stored in a .txt file by the 
#' `ps_write_bg_to_file()` function. In this case, only information about the 
#' background mean and standard deviation are stored for each matrix. These 
#' metrics will be later used for the computation of z-score by the pscan 
#' function. 
#' 
#' If a full background PSMatrixList is required, that includes all the 
#' background scores for each oligonucleotide hits, their position, strand and 
#' names, the user should store the output by the save() function and set 
#' the fullBG flag as TRUE. Note that this process will generate a very large 
#' file, that can reach several gigabytes in size. 
#' 
#' Note: Currently, it is not possible to use this function with the 
#' Arabidopsis thaliana organism due to the transcript nomenclature used. 
#' The function is designed to handle transcript names without version 
#' extensions (e.g., NM_30287 instead of NM_30287.1). 
#' However, in Arabidopsis thaliana, the format typically follows a 
#' gene.transcript pattern (e.g., AT1G01010.1). 
#' This causes the function to truncate the transcript identifier, effectively 
#' removing essential transcript-level information. As a result, the function 
#' may return incorrect results for this organism. This is valid for any 
#' organism for which this type of format is used.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @return 
#' A `PSMatrixList` object, containing each motif matrix from `pfms`, 
#' background-scored against the sequences in `x`. 
#' 
#' @seealso \code{\link{pscan_fullBG}}, \code{\link{ps_write_bg_to_file}}
#' 
#' @examples
#' # Note that the generation of the example may take few minutes
#' #
#' # Load the example dataset for promoter sequences (hg38 assembly, 
#' # -200 +50 bp in respect to the TSS).
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Load the example dataset for JASPAR2020 matrices collection for 
#' # vertebrates.
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#'
#' # Generate the background-scored motif matrices
#' bg_matrices <- ps_build_bg(prom_seq, J2020, 
#'                            BPPARAM = BiocParallel::SnowParam(1))
#' # Use BiocParallel::MulticoreParam() for Unix like systems. 
#' bg_matrices
#' bg_matrices[[1]]
#' 
#' # Example for full-background generation
#' full_bg_matrices <- ps_build_bg(prom_seq, J2020, 
#'                            BPPARAM = BiocParallel::SnowParam(1), 
#'                            fullBG = TRUE)
#'                            
#' full_bg_matrices[[1]]
#' 
#' @export
ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions(), 
                        fullBG = FALSE)
{
  .ps_checks(x, pfms, type = 1)
  
  x_unique <- BiocGenerics::unique(x)
  x_unique <- .clean_sequence(x_unique)
  
  pfms <- as(pfms, "PSMatrixList")

  pfms <- bplapply(
    pfms, 
    FUN = ps_scan, 
    x_unique, 
    BG = TRUE,
    fullBG = fullBG,
    BPPARAM=BPPARAM, 
    BPOPTIONS = BPOPTIONS
  )
  
  pfms <- do.call(PSMatrixList, pfms)
  
  if(fullBG == TRUE)
    pfms <- .mapping_unique_names(x, pfms)
  
  return(pfms)
  
}

#' Import background-scored matrices from a file 
#' 
#' Reconstructs a `PSMatrixList` object containing background-scored motif 
#' matrices, using a previously computed background dataset stored in a file. 
#' The input file must have been generated by the `ps_write_bg_to_file()` 
#' function from a background `PSMatrixList` object.
#' 
#' @param file A character string with the path to the input file. 
#'   This file contains the score distribution statistics for the 
#'   frequency matrices computed on a 
#'   background set of regulatory sequences.
#'   The file should be in tabular format, where the first column contains 
#'   the frequency matrix identifiers (row names) and subsequent columns 
#'   contain the background statistics 
#'   (e.g., \code{BG_SIZE}, \code{BG_MEAN}, \code{BG_STDEV}).
#'   Background statistics are generated with the `ps_build_bg` function and 
#'   can be written to a file with `ps_write_bg_to_file`.
#' 
#' @param pfms A `PFMatrixList` object containing position frequency matrices 
#' representing transcription factor binding preferences, obtained, 
#' for example, from the JASPAR database. 
#' The provided list must correspond to the same 
#' frequency matrix list used to build the background stored in `file`.
#'
#' @details 
#' This function:
#' \itemize{
#'    \item Validates the input type: must be a character string.  
#'    \item Reads the input file into a `data.frame`, with row names taken 
#'    from the first column. 
#'    \item Calls `ps_build_bg_from_table` with the created `data.frame` 
#'    and `pfms` as input.
#' }
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These file are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples. 
#' 
#' Other background datasets are aviable at the public repository
#' PscanR_background on GitHub: 
#' \url{https://github.com/dianabetelli/PscanR_backgrounds}
#' See vignettes for further details on the type of background available.
#' 
#' @seealso \code{\link{ps_build_bg}}, \code{\link{ps_write_bg_to_file}}, 
#' \code{\link{ps_build_bg_from_table}}
#' 
#' @return 
#' A `PSMatrixList` object, containing each motif matrix from `pfms`, 
#' background-scored using the values provided in `file`. 
#' See `ps_build_bg_from_table` for more details.
#' 
#' @examples
#' # Load a background information file
#' file_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                          package = "PscanR")
#'
#' # Load the example dataset for JASPAR2020 matrices collection for 
#' # vertebrates.
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#'
#' # Generate the background-scored motif matrices from file
#' bg_matrices <- ps_retrieve_bg_from_file(file_path, J2020)
#' bg_matrices
#' bg_matrices[[2]]
#' @seealso \code{\link{ps_build_bg_from_table}}
#' 
#' @export
ps_retrieve_bg_from_file <- function(file, pfms)
{
  .ps_checks(file, pfms, type = 2)
  
  short.matrix <- read.table(file, header = FALSE, row.names = 1, skip = 1)
  
  colnames(short.matrix)[seq_len(3)] <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
  
  pfms <- ps_build_bg_from_table(short.matrix, pfms)
}

#' Apply background parameters to position frequency matrices
#' 
#' Associates background parameters from a `data.frame` with a list of position 
#' frequency matrices (PFMs) and returns the updated matrices as 
#' a `PSMatrixList` object.
#' 
#' @param x A `data.frame` containing background parameters for each motif.
#'    The data frame should include the following columns:
#'    \itemize{
#'       \item `BG_SIZE`: Background set size.
#'       \item `BG_MEAN`: Mean score of background hits.
#'       \item `BG_STDEV`: Standard deviation of background hits score.
#'    }
#'    The row names of `x` should correspond to the identifiers of the motifs 
#'    in `pfms`.
#'    
#' @param pfms A `PFMatrixList` object containing position frequency matrices 
#'   representing transcription factor binding preferences, obtained, for 
#'   example, from the JASPAR database.
#'   The number of elements in `pfms` should match the number of row of `x`, 
#'   and their identifiers should correspond.
#' 
#' @details 
#' This function: 
#' \itemize{
#'    \item Validates the input type. `x` must be a `data.frame`.
#'    \item Converts each elements of pfms into `PSMatrix` objects. 
#'    \item A warning is issued if the number of elements in `pfms` doesn't 
#'    match the number of row of `x`.
#' }
#' 
#' Note: This function does not compute new background scores — it assigns 
#' existing background parameters to each motif, assuming they were previously 
#' computed.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @return 
#' A `PSMatrixList` object containing each motif from `pfms`, 
#' enriched with background scoring parameters from `x`.
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
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#' J2020_subset <- J2020[1:3] # match the number of rows in `backgound_data`
#' 
#' rownames(background_data) <- c("MA0004.1", "MA0006.1", "MA0019.1")
#' 
#' # Generate background-scored motif matrices
#' bg_matrices <- ps_build_bg_from_table(background_data, J2020_subset)
#' bg_matrices
#' bg_matrices[[1]]
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

#' Extract background statistics from a `PSMatrixList` object
#' 
#' Extracts background statistics (size, mean, and standard deviation) 
#' from a `PSMatrixList` of Position Weight Matrices representing 
#' transcription factor binding preferences, and returns a data frame 
#' containing these metrics.
#'
#' @param pfms A `PSMatrixList` of Position Weight Matrices representing 
#'   transcription factor binding preferences . 
#'   Each element must be a `PSMatrix` object or coercible to one.
#'
#' @return 
#' A `data.frame` with one row per matrix (corresponding to a 
#' transcription factor) in `pfms`. 
#' The row names will match the matrix identifiers.
#' 
#' Columns: 
#' \itemize{
#'   \item `BG_SIZE`: An integer vector containing the background set size for 
#'   each PSM. 
#'   \item `BG_MEAN`: A numeric vector containing the mean of the background 
#'   hits score for each PSM.
#'   \item `BG_STDEV`: A numeric vector containing the standard deviation 
#'   of the background hits score for each PSM. 
#' }
#' 
#' @details
#' 
#' This function is useful for exporting background-scored matrices to a 
#' tabular format. The resulting `data.frame` can be saved and later 
#' reapplied to PFMs using `ps_build_bg_from_table()`.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @seealso \code{\link{ps_build_bg_from_table}}, 
#' \code{\link{ps_retrieve_bg_from_file}}
#' 
#' @examples
#' # Retrieve motif matrices for vertebrates from JASPAR2020
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
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

#' Save background statistics from a `PSMatrixList` object to a file
#'
#' Saves background statistics (such as size, mean, and standard deviation) 
#' for each Position Weight Matrix in a `PSMatrixList` object to a specified 
#' file. 
#'
#' @param pfms A `PSMatrixList` object of Position Weight Matrices 
#'   representing transcription factor binding preferences, obtained, for 
#'   example, from the JASPAR database. 
#'   Each matrix should already include background statistics 
#'   (e.g., computed from all promoter regions in the organism of 
#'   study). Each element must be a `PSMatrix` object or coercible to one.
#'
#' @param file A character string specifying the path to the output file where 
#'   the background statistics should be saved.
#'
#' @details 
#' 
#' This function first extracts background statistics from the input `pfms` 
#' using `ps_get_bg_table()`. It then writes the result to the specified file 
#' in a tab-delimited format. The output file begins with a header line:
#' `[SHORT TFBS MATRIX]`. This header is used internally for compatibility with 
#' `ps_retrieve_bg_from_file()`, which expects it when reading.
#' 
#' The output does not include column names; the first column contains matrix 
#' identifiers, followed by `BG_SIZE`, `BG_MEAN`, and `BG_STDEV`.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples. 
#'
#' @return None. It saves the given background statistics to the specified file
#' in a tab-delimited format.
#'
#' @examples
#' # Since the function create a .txt file in the user working directory, 
#' # this example will not run automatically.
#' # User can test it in the console.
#' \dontrun{
#' # Load the example dataset for JASPAR2020 matrices collection 
#' # for vertebrates.
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
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
#' }
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

#' Generate PSMatrixList from Background data and JASPAR Matrix
#' 
#' This function retrieves background data for a specified organism and 
#' promoter region, and then fetches a corresponding matrix set from a 
#' specified JASPAR version. 
#' It returns a `PSMatrixList` object with background statistics.
#'
#' @param JASPAR_matrix A character string specifying the JASPAR database 
#'    version. You can choose between 'JASPAR2020', 'JASPAR2022', or 
#'    'JASPAR2024' (non case sensitive). 
#' @param org A string representing the organism acronym. Accepted values are: 
#'    `hs` (Homo sapiens),
#'    `mm` (Mus musculus), 
#'    `at` (Arabidopsis thaliana),
#'    `sc` (Saccharomyces cerevisiae),
#'    `dm` (Drosophila melanogaster).
#' @param prom_reg A numeric vector of two integers as `[upstream, downstream]`
#'    representing the promoter region relative to the TSS. 
#'    Allowed combinations:
#'    \itemize{
#'      \item 200 base pairs upstream and 50 downstream base pair from the TSS 
#'      (c(-200, 50))
#'      \item 450 bp upstream and 50 downstream (c(-450, 50))
#'      \item 500 bp upstream and 0 downstream (c(-500, 0))
#'      \item 950 bp upstream and 50 downstream (c(-950, 50))
#'      \item 1000 bp upstream and 0 downstream (c(-1000, 0))
#'      } 
#' @param assembly A string representing the assembly version for Human or 
#'   mouse. For `"hs"` you can choose between `"hg38"` or the latest `"hs1"`.
#'   For `"mm"` you can specify `"mm10"` or `"mm39"`. Default is character().
#' @param version A string indicating the version number of the desired 
#'   background. Since annotations may evolve over time, multiple versions 
#'   of the same background could exist. Default is '1', corresponding 
#'   to the first set of backgrounds generated between late 2024 and early 2025. 
#'   Currently, only this initial version is available. 
#' @param destfile A string indicating the path where the downloaded background 
#'   .txt file should be saved. This tab-separated file contains the matrix 
#'   identifiers, background size, average background score, and standard 
#'   deviation. See \code{\link{ps_retrieve_bg_from_file}} for details on how to 
#'   use this file. Default is NULL, so the file is not saved in user's working 
#'   environment. 
#'   
#' @details
#' The background files are downloaded from the GitHub repository:
#' \url{https://github.com/dianabetelli/PscanR_backgrounds}
#' This function automatically fetches the appropriate background file and 
#' combines it with the specified JASPAR matrix collection to create a 
#' `PSMatrixList` with background statistics.
#'
#' @return A `PSMatrixList` object with background-scored motif matrices.
#' 
#' @export
#'
#' @examples
#' # Note: when running this example, you may see a message indicating that 
#' # a file is being downloaded. This is expected behavior and not an error — 
#' # it simply informs you that background data is being retrieved.
#' bg_matrices <- generate_psmatrixlist_from_background('Jaspar2020', 'hs', 
#'                                                      c(-200,50), 'hg38')
#' bg_matrices
#' bg_matrices[[4]]
#' 
#' @import httr
#' @importFrom TFBSTools getMatrixSet
generate_psmatrixlist_from_background <- function(JASPAR_matrix, org, prom_reg, 
                                                  assembly = character(),
                                                  version = '1',
                                                  destfile = NULL){
  
  organism_map <- c("hs" = assembly, "mm" = assembly, "at" = "TAIR9", 
                       "sc" = "sacCer3", "dm" = "dm6")
  tax_map <- c('hs' = 'vertebrates', 'mm' = 'vertebrates', 'at' = 'plants',
                  'sc' = 'fungi','dm' = 'insects') 
  
  org_assembly <- if (org %in% names(organism_map)) organism_map[[org]] else NULL
  
  J_name <- toupper(JASPAR_matrix)

  Jversion <- substr(JASPAR_matrix,7, 10)
  
  p_up <- abs(prom_reg[1])
  p_down <- abs(prom_reg[2])
  
  file_suffix <- ifelse(org_assembly == "TAIR9", 
                        sprintf("TAIR.psbg%s.txt", version),
                        sprintf("UCSC.psbg%s.txt", version))
  
  file_name <- paste0('J', Jversion, '_', org_assembly, '_', p_up, 'u_', p_down, 
                      'd_', file_suffix)
  
  AvailableBG <- get_availableBG()
  
  if(!(file_name %in% AvailableBG))
    stop(sprintf("Invalid file name: '%s'\n\nPlease run 'get_availableBG()' to see the list of available background files", 
                 file_name))
  
  BG_path <- .download_background(file = file_name, destfile = destfile)
  
  opts <- list("collection" = "CORE", "tax_group" = tax_map[[org]])
  
  J_matrix <- switch(J_name,
    'JASPAR2020' = TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts),
    'JASPAR2022' = TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts),
    'JASPAR2024' = {
      httr::set_config(httr::config(ssl_verifypeer = 0L))
      JASPAR2024 <- JASPAR2024::JASPAR2024()
      JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), 
                                          JASPAR2024::db(JASPAR2024))
      J_matrix <- TFBSTools::getMatrixSet(JASPARConnect, opts)
    }
  )
  
  ps_retrieve_bg_from_file(BG_path, J_matrix)
}

#' Get Available Pre-Computed Background Files 
#' 
#' This function retrieves the names of all available .txt files containing 
#' background information, that serves for the construction of the background 
#' PSMatrixList. If a keyword is specified, the function will return only the 
#' file names that contains that string. Note that the filtering parameter 
#' works only if you know the nomenclature used to name the file. See the 
#' details paragraph for further informations. 
#' 
#' @param keyword A string to filter the file names. Default is NULL, so the 
#'    complete list of available file names is printed
#' 
#' @details
#' Some information for the filtering: 
#' \itemize{
#'    \item Type of collection: To filter by JASPAR version, 
#'    search by year prefixed by 'J', e.g., 'J2020'.
#'    \item To filter by genome assembly, use assembly identifiers like 
#'    'mm10' or 'mm39'.
#'    \item To filter by promoter region, use strings like '500u_0d', where 'u' 
#'    stands for 'upstream' and 'd' for 'downstream'
#'    \item You may also filter by version number, e.g., '1.txt'.
#'    \item Multiple keywords can be combined (e.g., 'hs1.*450u_50d').}
#' 
#' @seealso \code{\link{generate_psmatrixlist_from_background}}, \code{\link{ps_retrieve_bg_from_file}}
#' 
#' @examples
#' head(get_availableBG())
#' get_availableBG('mm10')
#' get_availableBG('hs1.*450u_50d') 
#' 
#' @import httr
#' @export
get_availableBG <- function(keyword = NULL) {
  
  url <- 'https://api.github.com/repos/dianabetelli/PscanR_backgrounds/contents/BG_files'
  
  response <- httr::GET(url)
  
  httr::stop_for_status(response)
  
  files_info <- content(response)
  file_names <- sapply(files_info, function(file) file$name)
  
  if(!is.null(keyword)){
    if(!is.character(keyword))
      stop('Keyword parameter must be a string')
    filtered_file_names <- file_names[grep(keyword, file_names)]
    if(length(filtered_file_names)!=0)
      return(filtered_file_names)
    else
      stop('Found 0 match with: ', keyword)
  }
  else
    return(file_names)
}