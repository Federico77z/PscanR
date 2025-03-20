#' @import TFBSTools 
#' @importFrom TFBSTools PFMatrix
.PSMatrix <- setClass("PSMatrix", 
                      slots = representation(ps_bg_avg="numeric",
                                             ps_fg_avg="numeric",
                                             ps_bg_std_dev="numeric", 
                                             ps_bg_size="integer",
                                             ps_fg_size="integer",
                                             ps_hits_pos="integer", 
                                             ps_hits_pos_bg = "integer",
                                             ps_hits_strand="character",
                                             ps_hits_strand_bg = "character",
                                             ps_hits_score="numeric",
                                             ps_hits_score_bg = "numeric",
                                             ps_hits_oligo="character",
                                             ps_hits_oligo_bg="character",
                                             ps_zscore="numeric",
                                             ps_pvalue="numeric",
                                             ps_seq_names="character",
                                             ps_bg_seq_names="character",
                                             .PS_PSEUDOCOUNT="numeric",
                                             .PS_ALPHABET="integer"),
                      contains="PFMatrix")

#' An S4 Class for a Single Positiont-Specific Matrix
#'
#' This function creates a `PSMatrix` object, which contains a single 
#' Position Frequency Matrix (PFM) and associated metadata related to promoter 
#' sequence motif hits, such as background and foreground averages, standard 
#' deviation, size, z-scores, and other related information. 
#'
#' @param pfm An object of class `PFMatrix` from the `TFBSTools` package, 
#'     typically a position frequency matrix (PFM) representing the motif for 
#'     which all the organism's promoter regions will be scanned.  
#' @param ps_bg_avg Numeric. The background average value, default = `NA`.
#'    It is the average binding score of promoter regions when scanned with a 
#'    PWM.
#' @param ps_fg_avg Numeric. The foreground average value, default = `NA`.
#'    Represents the average binding score of a set of promoter regions of 
#'    co-expressed or co-regulated genes when scanned with a PWM.
#' @param ps_bg_std_dev Numeric. The background standard deviation of PWM 
#'    scores, default = `NA`.
#' @param ps_bg_size Integer. The size of the background promoter region, 
#'    default = `NA`.
#' @param .PS_PSEUDOCOUNT Numeric. The Pseudo-Count to add to avoid division
#'    by zero. Default = `0.01`
#' @param ... Additional arguments passed to other methods or used in the 
#'    initialization.
#'    
#' @return A `PSMatrix` object containing the provided statistics, initialized 
#'    with default or provided values. 
#' 
#' @details
#' When a `PSMatrix` object is created or modified, it is automatically 
#' validated using an internal function. The validation ensures that the 
#' background and foreground statistics are properly formatted, and that motif 
#' hit-related vectors have consistent lengths.
#'    
#' @export
#'
#' @examples
#'
#'# Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#' 
#' result <- PSMatrix(
#'   pfm = J2020[[1]],
#'   ps_bg_avg = 0.25,        
#'   ps_fg_avg = 0.5,         
#'   ps_bg_std_dev = 0.05,    
#'   ps_bg_size = 250L        
#'   )
#' print(result)
PSMatrix <- function(pfm, ps_bg_avg = as.numeric(NA), 
                     ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
                     ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01, ...)
{
  #  .ps_required_packages()
  .ps_norm_matrix(.PSMatrix(pfm, ps_bg_avg = ps_bg_avg, 
                            ps_fg_avg = ps_fg_avg, 
                            ps_bg_std_dev = ps_bg_std_dev, 
                            ps_bg_size = ps_bg_size, 
                            ps_fg_size = as.integer(NA), 
                            ps_zscore = as.numeric(NA), 
                            ps_pvalue = as.numeric(NA), 
                            ps_seq_names = character(),
                            ps_bg_seq_names = character(),
                            .PS_PSEUDOCOUNT = .PS_PSEUDOCOUNT, 
                            ps_hits_pos = integer(),
                            ps_hits_pos_bg = integer(),
                            ps_hits_strand = character(),
                            ps_hits_strand_bg = character(),
                            ps_hits_score = numeric(),
                            ps_hits_score_bg = numeric(), 
                            ps_hits_oligo = character(),
                            ps_hits_oligo_bg = character(),
                            .PS_ALPHABET = setNames(seq_len(4), 
                                                    c("A","C","G","T"))))
}

#' @importFrom TFBSTools PFMatrixList
.PSMatrixList <-setClass("PSMatrixList", contains ="PFMatrixList", 
                         slots = list(transcriptIDLegend = 'character'))

#' An S4 Class for Storing Position-Specific Matrices
#' 
#' This function creates a `PSMatrixList` object, which is a container for 
#' managing multiple `PSMatrix` object. It is similar to `PFMatrixList` but 
#' extends its functionality.
#'
#' @param ... Objects of class `PSMatrix` to include in the list.
#' @param use.names Logical. Assert whether to use names from the input objects.
#'    Default = `TRUE`
#' @param transcriptIDLegend Named character vector. The names correspond to the
#'   IDs of all the transcript expressed in the organism of study, whereas 
#'   corresponding values are the IDs of transcript retained by the unique() 
#'   function when multiple identical sequences exist. So, if multiple identical 
#'   sequences exist (e.g., ID1, ID2, ID3, and ID4), and unique() retains only 
#'   ID2, the mapping will associate each original name with its unique 
#'   counterpart (ID1 → ID2, ID2 → ID2, ID3 → ID2, ID4 → ID2). 
#'
#' @return A `PSMatrixList` object, which is a list containing `PSMatrix` 
#'    objects. Each element in the list corresponds to a `PSMatrix` object 
#'    provided as input.
#'
#' @examples
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
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
#' result <- PSMatrixList(PSM1, PSM2)
#' ps_results_table(result)
#' @export
PSMatrixList <- function(..., transcriptIDLegend = character(), use.names = TRUE)
{
  listData <- list(...)
  XMatrixList(listData, 
              use.names = use.names, 
              type = "PSMatrixList", 
              matrixClass = "PSMatrix")
}