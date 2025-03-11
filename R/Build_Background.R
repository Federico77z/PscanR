#' Build background matrices in parallel
#' 
#' Generates a background probability profile matrix list using the BiocParallel
#' framework for parallel execution. It computes background scores for each 
#' motif matrix, collected, for example, from the JASPAR database, by scanning 
#' them against a set of regulatory sequences (e.g., gene promoters) using the 
#' `ps_scan` function.
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
#' @param megaBG Logical. Default is FALSE. When set to TRUE, it creates a 
#'   mapping between all sequence names in the organism of study 
#'   and the corresponding names retained after applying the unique() function. 
#'   For example, if multiple identical sequences exist (e.g., ID1, ID2, ID3, 
#'   and ID4), and unique() retains only ID2, the mapping will associate each 
#'   original name with its unique counterpart (ID1 → ID2, ID2 → ID2, ID3 → ID2, 
#'   ID4 → ID2). This vector is stored in the transcriptIDLegend slot of the 
#'   PSMatrixList and helps generate a complete background PSMatrixList, 
#'   improving computational efficiency for the pscan_bg() function.
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
#' If a full background PSMatrixList is required, that includes all the 
#' background scores for each oligonucleotide hits, their position, strand and 
#' names, the user should store the output by the save() function and set 
#' the megaBG flag as TRUE. Note that this process will generate a very large 
#' file, that can reach several gigabytes in size. 
#' 
#' @return 
#' A `PSMatrixList` object, containing each motif matrix from `pfms`, 
#' background-scored against the sequences in `x`. 
#' 
#' @seealso \code{\link{pscan_bg}}, \code{\link{ps_write_bg_to_file}}
#'
#' @importFrom txdbmaker makeTxDbFromUCSC 
#' @importFrom GenomicFeatures promoters
#' @importFrom Biostrings getSeq
#' @importFrom TFBSTools getMatrixSet
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import JASPAR2020
#' 
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#'
#' # Generate the background-scored motif matrices
#' bg_matrices <- ps_build_bg(prom_seq, J2020, 
#'                            BPPARAM = BiocParallel::SnowParam(1))
#' # Use BiocParallel::MulticoreParam() for Unix like systems. 
#' bg_matrices
#' bg_matrices[[1]]
#' 
#' # Example for mega-background generation
#' mega_bg_matrices <- ps_build_bg(prom_seq, J2020, 
#'                            BPPARAM = BiocParallel::SnowParam(1), megaBG = TRUE)
#'                            
#' mega_bg_matrices[[1]]
#' 
#' @export
ps_build_bg <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions(), 
                        megaBG = FALSE)
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
    megaBG = megaBG,
    BPPARAM=BPPARAM, 
    BPOPTIONS = BPOPTIONS
  )
  
  pfms <- do.call(PSMatrixList, pfms)
  
  if(megaBG == TRUE)
    pfms <- .mapping_unique_names(x, pfms)
  
  return(pfms)
  
}

#' Import background statistics from a file 
#' 
#' Reads and imports background score distribution statistics from a file, 
#' corresponding to a PSMatrixList object previously computed on a background 
#' set of regulatory sequences (e.g., the promoters of all transcripts in the 
#' organism of study, defined by coordinates relative to the TSS) by the 
#' `ps_build_bg()` function. `ps_retrieve_bg_from_file()` reads the background 
#' statistics stored into a file by `ps_write_bg_to_file()`. 
#' 
#' @param file A character string with the path to the input file. 
#'   This file contains the score distribution statistics for the 
#'   frequency matrices computed on a 
#'   background set of regulatory sequences.
#'   The file should be in tabular format, where the first column contains 
#'   the frequency matrix identifiers (row names) and subsequent columns contain 
#'   the background statistics (e.g., `BG_SIZE`, `BG_MEAN`, `BG_STDEV`).
#'   Background statistics are generated with the `ps_build_bg` function and can 
#'   be written to a file with `ps_write_bg_to_file`.
#' 
#' @param pfms A `PFMatrixList` object containing position frequency matrices 
#' representing transcription factor binding preferences, obtained, for example, 
#' from the JASPAR database. The provided list must correspond to the same 
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
#' @seealso \code{\link{ps_build_bg}}, \code{\link{ps_write_bg_to_file}}, 
#' \code{\link{ps_build_bg_from_table}}
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
#' file_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                          package = "PscanR")
#'
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
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
#' Associates background parameters from a data frame with a list of position 
#' frequency matrices (PFMs) and returns the list of updated matrices as 
#' `PSMatrixList` object.
#' 
#' @param x A `data.frame` containing background parameters for each motif.
#'    The data frame should have the following columns:
#'    \itemize{
#'       \item `BG_SIZE`: Background size.
#'       \item `BG_MEAN`: Mean of the background hits score.
#'       \item `BG_STDEV`: Standard deviation of the background hits score.
#'    }
#'    The number of rows in `x` must match the number of elements in `pfms`.
#'    
#' @param pfms A `PFMatrixList` object of position frequency matrices 
#'   representing transcription factor binding preferences, obtained, for 
#'   example, from the JASPAR database.
#'   The number of elements in `pfms` should match the number of row of `x`.
#' 
#' @details 
#' This function: 
#' \itemize{
#'    \item Validates the input type. `x` must be a `data.frame`.
#'    \item Converts each elements of pfms into `PSMatrix` class. 
#'    \item A warning is issued if the number of elements in `pfms` doesn't 
#'    match the number of row of `x`.
#' }
#' 
#' @return 
#' A `PSMatrixList` object containing each motif matrix from `pfms`, 
#' scored with background parameters provided by `x`.
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
#' transcription factor binding preferences and generates a table containing 
#' these metrics.
#'
#' @param pfms A `PSMatrixList` of Position Weight Matrices representing 
#'   transcription factor binding preferences . 
#'   Each element should be a `PSMatrix` object 
#'   (or should be coercible to `PSMatrix`).
#'
#' @return 
#' A `data.frame` with one row for each PFM in `pfms`.
#' 
#' Columns: 
#' \itemize{
#'   \item `BG_SIZE`: An integer vector representing the background size for 
#'   each PFMs. 
#'   \item `BG_MEAN`: A numeric vector representing the mean of the background 
#'   hits score for each PFMs.
#'   \item `BG_STDEV`: A numeric vector representing the standard deviation 
#'   of the background hits score for each PFMs. 
#' }
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
#'   representing transcription factor binding preferences, obtained for example 
#'   from the JASPAR database, containing background statistics computed by 
#'   scanning each `PSMatrix` against the set of promoter regions of all 
#'   transcripts in the organism of study. The promoter regions analyzed may 
#'   vary. Each element should be a `PSMatrix` object or coercible to `PSMatrix`.
#'
#' @param file A character string specifying the path to the output file where 
#'   the background statistics should be saved.
#'
#' @details 
#' This function retrieves the background statistics (size, mean, and 
#' standard deviation) using `ps_get_bg_table()` from the input `PSMatrixList` 
#' object, after having validated the inputs. Then, it writes the result to a 
#' specified file. 
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
#' # ps_write_bg_to_file(PSMatrixList_J2020, file_path)
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
#' This function retrieves background data for a specified organism and promoter 
#' region, and then fetches a corresponding matrix set from a specified JASPAR 
#' version. It returns a `PSMatrixList` object.
#'
#' @param JASPAR_matrix A character string specifying the JASPAR database 
#'    version. You can choose between 'JASPAR2020', 'JASPAR2022', or 
#'    'JASPAR2024' (non case sensitive). 
#' @param org A string representing the organism acronym. Accepted values are: 
#'    `hs` (Homo sapiens), `mm` (Mus musculus), `at` (Arabidopsis thaliana), 
#'    `sc` (Saccharomyces cerevisiae), `dm` (Drosophila melanogaster). 
#' @param prom_reg A numeric vector of two integers as `[upstream, downstream]`
#'    bp from the transcription start site (TSS) indicating the range 
#'    of the promoter region. You can choose between:
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
#'   For `"mm"` you can specify `"mm10"` or `"mm39"`.
#'
#' @return A `PSMatrixList` object created from the specified background 
#' file and JASPAR matrix.
#' 
#' @export
#'
#' @examples
#' bg_matrices <- generate_psmatrixlist_from_background('Jaspar2020', 'hs', 
#'                                                      c(-200,50), 'hg38')
#' bg_matrices
#' bg_matrices[[4]]
#' 
generate_psmatrixlist_from_background <- function(JASPAR_matrix, org, prom_reg, 
                                                  assembly){
  
  organism_map <- c("hs" = assembly, "mm" = assembly, "at" = "TAIR9", 
                       "sc" = "sacCer3", "dm" = "dm6")
  tax_map <- c('hs' = 'vertebrates', 'mm' = 'vertebrates', 'at' = 'plants',
                  'sc' = 'fungi','dm' = 'insects') 
  
  org_assembly <- organism_map[[org]]
  if(is.null(org_assembly))
    stop("Invalid organism acronym. Choose between: 'hs', 'mm', 'at', 'sc', 
         'dm'.")
  
  valid_JASPAR <- c("JASPAR2020", "JASPAR2022", "JASPAR2024")
  J_name <- toupper(JASPAR_matrix)
  if(!(J_name %in% valid_JASPAR))
    stop("Invalid database. Choose between: 'JASPAR2020', 'JASPAR2022', 
    'JASPAR2024' (non case sensitive).")
  
  version <- substr(JASPAR_matrix,7, 10)
  
  if(length(prom_reg) != 2 || !is.numeric(prom_reg))
    stop("Promoter region must be a numeric vector with two integer.")
  p_up <- abs(prom_reg[1])
  p_down <- abs(prom_reg[2])
  
  file_suffix <- ifelse(org_assembly == "TAIR9", "TAIR.psbg.txt", 
                        "UCSC.psbg.txt")
  
  file_name <- paste0('J', version, '_', org_assembly, '_', p_up, 'u_', p_down, 
                      'd_', file_suffix)
  
  BG_path <- system.file("extdata/BG_scripts", file_name, package = 'PscanR')
  
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
