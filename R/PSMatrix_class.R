#' An S4 class to represent a `PFMatrix` object with computed foreground and
#' background statistics 
#' 
#' The `PSMatrix` class extends the `PFMatrix` class from the `TFBSTools`
#' package, adding slots for foreground and background statistics, Z-score, 
#' P-value, and more.
#' 
#' @slot ps_bg_avg Numeric. The background average value.
#' @slot ps_fg_avg Numeric. The foreground average value.
#' @slot ps_bg_std_dev Numeric. The standard deviation of the background.
#' @slot ps_bg_size Integer. The size of the background dataset. 
#' @slot ps_fg_size Integer. The size of the foreground dataset.
#' @slot ps_hits_pos Integer. Positions of hits. 
#' @slot ps_hits_strand Character. Strand information for hits. 
#' @slot ps_hits_score Numeric. Score information of hits. 
#' @slot ps_hits_oligo Character. Oligonucleotide sequence of hits.
#' @slot ps_zscore Numeric. The Z-score value.
#' @slot ps_pvalue Numeric. The P-Value.
#' @slot ps_seq_names Character. The sequence name.  
#' @slot .PS_PSEUDOCOUNT Numeric. The pseudocount value used in calculations. 
#' @slot .PS_ALPHABET Integer. From `1` to `4` to represent the DNA alphabet. 
#' 
#' @exportClass PSMatrix
#' @import TFBSTools 
#' @importFrom TFBSTools PFMatrix
.PSMatrix <- setClass("PSMatrix", 
                      slots = representation(ps_bg_avg="numeric",
                                             ps_fg_avg="numeric",
                                             ps_bg_std_dev="numeric", 
                                             ps_bg_size="integer",
                                             ps_fg_size="integer",
                                             ps_hits_pos="integer", 
                                             ps_hits_strand="character", 
                                             ps_hits_score="numeric",
                                             ps_hits_oligo="character",
                                             ps_zscore="numeric",
                                             ps_pvalue="numeric",
                                             ps_seq_names="character",
                                             .PS_PSEUDOCOUNT="numeric",
                                             .PS_ALPHABET="integer"), 
                      contains="PFMatrix")

#' Create a PSMatrix object
#' 
#' This function creates an instance of the `PSMatrix` class by initializing its
#' slots with default values. 
#'
#' @param pfm An object of class `PFMatrix` from the `TFBSTools` package 
#'    representing a position frequency matrix.  
#' @param ps_bg_avg Numeric. The background average value, default = `NA`.
#' @param ps_fg_avg Numeric. The foreground average value, default = `NA`.
#' @param ps_bg_std_dev Numeric. The standard deviation of the background, 
#'    default = `NA`.
#' @param ps_bg_size Integer. The size of the background dataset, 
#'    default = `NA`.
#' @param .PS_PSEUDOCOUNT Numeric. The pseudocount to add to avoid division
#'    by zero. Default = `0.01`
#' @param ... Additional arguments to pass to downstream function.
#'
#' @return An object of `PSMatrix` class.
#' @export
#'
#' @examples
#'
#'# Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
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
                            .PS_PSEUDOCOUNT = .PS_PSEUDOCOUNT, 
                            ps_hits_pos = integer(), 
                            ps_hits_strand = character(), 
                            ps_hits_score = numeric(), 
                            ps_hits_oligo = character(), 
                            .PS_ALPHABET = setNames(seq_len(4), 
                                                    c("A","C","G","T"))))
}

#' An S4 class to represent a `PSMatrixList` object 
#' 
#' The `PSMatrixList` class extend the `PFMatrixList` class and represents a
#' collection of `PSMatrix` object. It allows to aggregate and iterate on 
#' multiple `PSMatrix` as a cohesive group.
#' 
#' @details
#' This class extends the `PFMatrixList` class without adding new slots. 
#' See `PFMatrixList` class documentation.
#'
#' @exportClass PFMatrixList
#' @importFrom TFBSTools PFMatrixList
.PSMatrixList <-setClass("PSMatrixList", contains ="PFMatrixList")

#' Create a `PSMatrixList` object
#' 
#' This function creates a `PSMatrixList` object, which is a container for 
#' managing multiple `PSMatrix` object.
#'
#' @param ... Objects of class `PSMatrix` to include in the list.
#' @param use.names Logical. Assert whether to use names from the input objects.
#'    Default = `TRUE`
#'
#' @return An object of class `PSMatrixList`.
#'
#' @examples
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
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
#' result
#' @export
PSMatrixList <- function(..., use.names = TRUE)
{
  listData <- list(...)
  XMatrixList(listData, 
              use.names = use.names, 
              type = "PSMatrixList", 
              matrixClass = "PSMatrix")
}

#' Retrieve or Manipulate Properties of PSMatrix Class
#'
#' These generics provide a standardized interface for accessing or manipulating 
#' various properties of `PSMatrix` objects.
#'
#' @param x An object of class `PSMatrix`.
#' @param ... Additional arguments passed to specific methods.
#' 
#' @name ps_generics
NULL


#' @describeIn ps_generics 
#' Retrieve the Z-score from an object
#'
#' @return For `ps_zscore`: a numeric value representing the Z-score.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_zscore(pfm1)
#' 
#' @export
setGeneric("ps_zscore", function(x, ...) standardGeneric("ps_zscore"))

#' @describeIn ps_generics 
#' Retrieve the P-Value from an object
#'
#' @return For `ps_pvalue`: a numeric value representing the P-Value.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_pvalue(pfm1)
#' 
#' @export
setGeneric("ps_pvalue", function(x, ...) standardGeneric("ps_pvalue"))

#' @describeIn ps_generics 
#' Retrieve the average background value from an object 
#'
#' @return For `ps_bg_avg`: a numeric value representing the average background.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_avg(pfm1)
#' 
#' @export
setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' @describeIn ps_generics 
#' Retrieve the average foreground value from an object
#'
#' @return For `ps_fg_avg`: a numeric value representing average foreground.
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_fg_avg(pfm1)
#' 
#' @export
setGeneric("ps_fg_avg", function(x, ...) standardGeneric("ps_fg_avg"))

#' @describeIn ps_generics 
#' Retrieve the background standard deviation from an object
#'
#' @return For `ps_bg_std_dev`: a numeric value representing the background 
#'    standard deviation.
#'    
#' @examples 
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_std_dev(pfm1)
#' 
#' @export
setGeneric("ps_bg_std_dev", function(x, ...) standardGeneric("ps_bg_std_dev"))

#' @describeIn ps_generics 
#' Retrieve the background size from an object
#'
#' @return For `ps_bg_size`: an integer value representing the size of the 
#'    background dataset.
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_size(pfm1)
#' 
#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

#' @describeIn ps_generics 
#' Retrieve the foreground size from an object
#'
#' @return For `ps_fg_size`: an integer value representing the size pf 
#'    the foreground dataset.
#'    
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_fg_size(pfm1)
#' 
#' @export
setGeneric("ps_fg_size", function(x, ...) standardGeneric("ps_fg_size"))

#' @describeIn ps_generics 
#' Retrieve the hits size from an object
#'
#' @return For `ps_hits_size`: an integer value representing the hits size.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_size(pfm1)
#' 
#' @export
setGeneric("ps_hits_size", function(x, ...) standardGeneric("ps_hits_size"))

#' @describeIn ps_generics 
#' Retrieve the hits score from an object
#'
#' @return For `ps_hits_score`: a numeric vector representing the score for 
#'    each hit.
#' 
#' @examples 
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_score(pfm1)
#' 
#' @export
setGeneric("ps_hits_score", function(x, ...) standardGeneric("ps_hits_score"))

#' @describeIn ps_generics 
#' Retrieve the hits z scores from an object. 
#'
#' @return For `ps_hits_z`: a numeric vector containing the z scores.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_z(pfm1)
#' 
#' @export
setGeneric("ps_hits_z", function(x, ...) standardGeneric("ps_hits_z"))

#' @describeIn ps_generics 
#' Retrieve the strand of hits from an object
#'
#' @return For `ps_hits_strand`: a character vector representing the strand 
#'    of each hit.
#'    
#' @examples 
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_strand(pfm1)
#' 
#' @export
setGeneric("ps_hits_strand", function(x, ...) standardGeneric("ps_hits_strand"))

#' @describeIn ps_generics 
#' Retrieve the position of hits from an object
#'
#' @return For `ps_hits_pos`: an integer vector representing the position of
#'    each hit.
#'    
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_pos(pfm1)
#' 
#' @export
setGeneric("ps_hits_pos", function(x, ...) standardGeneric("ps_hits_pos"))

#' @describeIn ps_generics 
#' Retrieve oligonucleotide sequence of hits from an object
#'
#' @return For `ps_hits_oligo`: a character vector representing the 
#'    oligonucleotide sequence of each hit. 
#' 
#' @export
setGeneric("ps_hits_oligo", function(x, ...) standardGeneric("ps_hits_oligo"))

#' @describeIn ps_generics 
#' Retrieve the hits table from an object
#'
#' @return For `ps_hits_table`: a `data.frame` of hits.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_table(pfm1)
#' 
#' @export
setGeneric("ps_hits_table", function(x, ...) standardGeneric("ps_hits_table"))
# not sure

#' @describeIn ps_generics 
#' Retrieve the sequence name from an object
#'
#' @return For `ps_seq_names`: a character vector of sequence names.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_seq_names(pfm1)
#' 
#' @export
setGeneric("ps_seq_names", function(x, ...) standardGeneric("ps_seq_names"))

#' @describeIn ps_generics
#' Retrieve the pseudocount value
#' 
#' @return For `.PS_PSEUDOCOUNT`: a numeric value representing the pseudocount.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .PS_PSEUDOCOUNT(pfm1)
#' 
#' @export
setGeneric(".PS_PSEUDOCOUNT", function(x, ...) standardGeneric(".PS_PSEUDOCOUNT"))

#' @describeIn ps_generics
#' Retrieve the numbers representing the DNA alphabet
#'  
#' @return For `.PS_ALPHABET`: integer values.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .PS_ALPHABET(pfm1)
#' 
#' @export
setGeneric(".PS_ALPHABET", function(x, ...) standardGeneric(".PS_ALPHABET"))

#' Normalize an Object 
#' 
#' @param x `PSMatrix` Object 
#' @param ... Additional parameters
#' 
#' @return A normalized object. The specific class depends on the implementation
#'    of the method. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .ps_norm_matrix(pfm1)
#' 
#' @export
setGeneric(".ps_norm_matrix", function(x, ...) standardGeneric(".ps_norm_matrix"))

#' Retrieve the Sequences Names from an Object
#' 
#' A generic function designed to process or retrieve sequence names 
#' from an object.
#' 
#' @param x The object from which sequence names will be processed or retrieved.
#' @param out A placeholder 
#' 
#' @return Depends on the specific method implementation 
#' 
#' @export
setGeneric(".ps_seq_names", function(x, out) standardGeneric(".ps_seq_names"))

#' @describeIn ps_generics 
#' Perform a scan operation on an Object 
#'
#' @return A `data.frame` of hits. 
#' 
#' @export
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))

#' Retrieve values from a table 
#' 
#' @param x A `data.frame`
#' @param ... Additional parameters
#' 
#' @return A `data.frame` of values depending on the specific implementation 
#'    method.
#' 
#' @export
setGeneric(".ps_bg_from_table", function(x, ...) standardGeneric(".ps_bg_from_table"))

#' Generic Function for Scanning 
#' 
#' A generic function that performs scanning operations on the provided object.
#' 
#' @param x The object on which the scanning operation is performed.
#' @param ... Additional arguments passed to methods for `.ps_scan_s`.
#' 
#' @return The return value depends on the specific method implementation.
#' 
#' @export
setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))

#' Generic Function for Normalizing Scores
#' 
#' A generic function that computes normalized scores for the provided object.
#' 
#' @param x An object on which the normalization of scores is performed.
#' @param ... Additional arguments passed to methods for `.ps_norm_score`.
#' 
#' @return A numeric vector containing the normalized scores.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .ps_norm_score(pfm1)
#' 
#' @export
setGeneric(".ps_norm_score", function(x, ...) standardGeneric(".ps_norm_score"))

#setGeneric(".ps_assign_score", function(x, ...) standardGeneric(".ps_assign_score"))

#' Generic Function for Adding Hits
#' 
#' A generic function designed to add hits (matching elements) to an object.
#' 
#' @param x An object on which the hits are added.
#' @param ... Additional arguments passed to methods for `.ps_add_hits`.
#' 
#' @return An extended object. 
#' 
#' @export
setGeneric(".ps_add_hits", function(x, ...) standardGeneric(".ps_add_hits"))

#' Get Background Average Score
#'
#' Retrieves the background average score stored in a `PSMatrix` object.
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the background average score.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_avg(pfm1) 
#'
#' @export
setMethod("ps_bg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_avg
  
  return(out)
})

#' Get Foreground Average Score
#'
#' Retrieves the foreground average score stored in a `PSMatrix` object.
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the foreground average score. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_fg_avg(pfm1)
#' 
#' @export
setMethod("ps_fg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_avg
  
  return(out)
})

#' Get Z-Score
#' 
#' Retrieves the Z-Score stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in 
#'    the output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the `PSMatrix` Z-score.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_zscore(pfm1)
#'
#' @export
setMethod("ps_zscore", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_zscore
  
  return(out)
})

#' Get P-Value
#' 
#' Retrieves the P-Value stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the `PSMatrix` P-Value.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_pvalue(pfm1)
#' 
#' @export
setMethod("ps_pvalue", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_pvalue
  
  return(out)
})

#' Get The Sequences of Hits
#' 
#' Retrieves the oligonucleotide sequences of hits.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A character containing the sequences of hits of a `PSMatrix` object. 
#' 
#' @export
setMethod("ps_hits_oligo", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_oligo

  out <- .ps_seq_names(out, x)
  
  return(out)
})

#' Method for Extracting Sequence Names
#' 
#' This method extract and set sequences name for object of class `PSMatrix`.
#' 
#' @param x An object of class `PSMatrix`. 
#' @param out The object to which sequence names are assigned.
#' 
#' @return An object of class `PSMatrix` to which sequence names are assigned. 
#' 
#' @export
setMethod(".ps_seq_names", "PSMatrix", function(x, out) {
  
  if(!any(is.na(x@ps_seq_names)))
    names(out) <- x@ps_seq_names
  
  return(out)
})

#' Get Background Standard Deviation
#' 
#' Retrieves the background standard deviation stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric value representing the `PSMatrix` background standard 
#'    deviation.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_std_dev(pfm1)
#' 
#' @export
setMethod("ps_bg_std_dev", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_std_dev
  
  return(out)
})

#' Get Background Size
#' 
#' Retrieves the background dimension stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer value representing the `PSMatrix` background size.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_bg_size(pfm1)
#' 
#' @export
setMethod("ps_bg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_size
  
  return(out)
})

#' Get Foreground Size
#' 
#' Retrieves the foreground dimension stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer value representing the `PSMatrix` foreground size.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_fg_size(pfm1)
#' 
#' @export
setMethod("ps_fg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_size
  
  return(out)
})

#' Compute Hits Size
#' 
#' Calculates the hits dimension stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer value representing the `PSMatrix` hits size.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_size(pfm1)
#' 
#' @export
setMethod("ps_hits_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- length(x@ps_hits_pos)
  
  return(out)
})

#' Get Hits Score
#' 
#' Retrieves the hits scores stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric vector representing the `PSMatrix` hits scores.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_score(pfm1)
#' 
#' @export
setMethod("ps_hits_score", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_score
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Compute Z-Scores for Hits in a PSMatrix Object
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric vector of z-score of hits.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_z(pfm1)
#' 
#' @export
setMethod("ps_hits_z", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- (x@ps_hits_score - x@ps_bg_avg) / x@ps_bg_std_dev
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get the hits strand in a PSMatrix Object
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A character vector of `-` and `+` indicating the hits strand.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_strand(pfm1)
#' 
#' @export
setMethod("ps_hits_strand", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_strand
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})


#' Get Adjusted Positions for Hits in a `PSMatrix` Object
#' 
#' Retrieves the positions of hits stored in a `PSMatrix` object, optionally 
#' applying a position shift.
#' 
#' @param x A `PSMatrix` object.
#' @param pos_shift Integer. Specifies the amount to shift the position. 
#'    Default is set to `0`.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer vector containing the position of hits of a `PSMatrix` object.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_pos(pfm1)
#' 
#' @export
setMethod("ps_hits_pos", "PSMatrix", function(x, pos_shift = 0L, 
                                              withDimnames = TRUE) {
  out <- x@ps_hits_pos + pos_shift
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get the sequences name in a PSMatrix Object
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A character vector of names.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_seq_names(pfm1)
#' 
#' @export
setMethod("ps_seq_names", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_seq_names
  
  return(out)
})

#' Get the Pseudocount Value in a PSMatrix Object
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric values representing the pseudocount added for the calculations.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .PS_PSEUDOCOUNT(pfm1)
#' 
#' @export
setMethod(".PS_PSEUDOCOUNT", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_PSEUDOCOUNT
  
  return(out)
})

#' Get the Alphabet from a PSMatrix Object
#' 
#' This method retrieves the alphabet of nucleotides associated with a 
#' `PSMatrix` object. 
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A named vector (integer) representing the alphabet. The names of the 
#' vector are the nucleotide symbols (e.g., `A`, `C`, `G`, `T`), and the 
#' values are the corresponding indices (e.g., 1, 2, 3, 4).
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .PS_ALPHABET(pfm1)
#' 
#' @export
setMethod(".PS_ALPHABET", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_ALPHABET
  
  return(out)
})

#' Retrieve a PSHits Table with Scores, Positions, Strands, and Oligonucleotides
#' 
#' @param x An object of class `PSMatrix`. Should contain the following 
#'    slots:
#'    \itemize{
#'      \item `ps_hits_score`: a numeric vector of hit scores.
#'      \item `ps_hits_strand`: a character vector of strand information
#'      (`-` and `+`).
#'      \item `ps_hits_oligo`: a character vector of oligo sequences.}
#' @param pos_shift Integer. Value for which the position gets shifted.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'  
#' @return A `data.frame` ordered by decreasing score value, 
#' with the following columns:
#' \itemize{
#'   \item `SCORE`: the score value of each hit
#'   \item `POS`: the position value of each hit
#'   \item `STRAND`: the strand information for each hit
#'   \item `OLIGO`: the oligo sequence corresponding to each hit} 
#' Row names correspond to the sequences name. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' ps_hits_table(pfm1)
#' 
#' @export
setMethod("ps_hits_table", "PSMatrix", function(x, pos_shift = 0L, 
                                                withDimnames = TRUE) {
  
  out <- data.frame("SCORE" = x@ps_hits_score, 
                    "POS" = ps_hits_pos(x, pos_shift = pos_shift), 
                    "STRAND" = x@ps_hits_strand,
                    "OLIGO" = DNAStringSet(x@ps_hits_oligo), 
                    row.names = x@ps_seq_names)
  
  out <- out[with(out, order(SCORE, POS, decreasing = c(TRUE,FALSE))),]
  
  return(out)
})

#' Add Hits Information to a `PSMatrix` Object
#' 
#' The `.ps_add_hits` method is used to add hit information to a `PSMatrix` 
#' object, including positional, strand, and score data. 
#' 
#' @param x A `PSMatrix` object to which the hits information will be added.
#' @param Pos An integer vector representing the positions of the hits.
#' @param Strand A character vector indicating the strand 
#'    of the hits (e.g., "+" or "-").
#' @param Score A numeric vector containing the scores of the hits.
#' @param Oligo A character vector representing oligonucleotide sequences 
#'    associated with the hits.
#' @param BG A logical value indicating whether to treat the 
#'    hits as background (default: `FALSE`).
#' @param withDimnames  Logical, whether to include dimension names in 
#'    the output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return The updated `PSMatrix` object.   
#' 
#' @export
setMethod(".ps_add_hits", "PSMatrix", 
          function(x, Pos, Strand, Score, Oligo, BG = FALSE, 
                   withDimnames = TRUE) {
  
  x@ps_hits_pos <- Pos
  x@ps_hits_strand <- Strand
  x@ps_hits_score <- Score
  x@ps_hits_score <- .ps_norm_score(x)
  
  if(BG)
  {
    ps_bg_size(x) <- length(x@ps_hits_pos)
    ps_bg_avg(x) <- mean(x@ps_hits_score, na.rm = TRUE)
    ps_bg_std_dev(x) <- sd(x@ps_hits_score, na.rm = TRUE)
    
    x@ps_hits_pos <- integer()
    x@ps_hits_strand <- character()
    x@ps_hits_score <- numeric()
  }
  else
  {
    if(!is.na(x@ps_bg_avg) && !is.na(x@ps_bg_std_dev))
    {
      ztest <- z.test(x@ps_hits_score, 
                      mu = x@ps_bg_avg, 
                      sigma.x = x@ps_bg_std_dev, 
                      alternative = "greater")
      
      x@ps_zscore <- ztest$statistic["z"]
      x@ps_pvalue <- as.numeric(ztest$p.value)
      x@ps_fg_avg <- mean(x@ps_hits_score, na.rm = TRUE)
      x@ps_fg_size <- length(x@ps_hits_pos)
      x@ps_hits_oligo <- Oligo
    }
  }
  
  return(x)
})

#' Normalize the Hit Scores in a `PSMatrix` Object
#' 
#' This method computes normalized hit scores in a `PSMatrix` Object
#' 
#' @param x An object of class `PSMatrix`.
#' 
#' @return A numeric vector of normalized scores.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .ps_norm_score(pfm1)
#' 
#' @export
setMethod(".ps_norm_score", "PSMatrix", function(x) {
  
  ps_score <- 1 + ((maxScore(Matrix(x)) - ps_hits_score(x)) / 
         (minScore(Matrix(x)) - maxScore(Matrix(x))))
  
  return(ps_score)
})

#' Populates Background Statistics in a `PSMatrix` Object from a Table
#' 
#' This method updates a `PSMatrix` object with background statistics 
#' retrieved from a provided table.
#' 
#' @param x An object of class `PSMatrix`.
#' @param short.matrix A data.frame or matrix containing background statistics.
#'    The table must have row names matching `ID(x)` and columns `BG_SIZE`,
#'    `BG_MEAN`, `BG_STDEV`.
#' 
#' @details
#' This method uses the `ID(x)` as a key to locate the matching rows in the
#' `short.matrix` table and updates the following background fields in the 
#' `PSMatrix` object:
#' \itemize{
#'   \item `ps_bg_size`: the size of the background data set.
#'   \item `ps_bg_avg`: the mean value of the background scores.
#'   \item `ps_bg_std_dev`: the standard deviation of the background scores.}
#'   
#' If no matching row is found in the table, a warning is issued. 
#' 
#' @return The updated `PSMatrix` object.
#' 
#' @export
setMethod(".ps_bg_from_table", "PSMatrix", function(x, short.matrix) {
  
  if(any(row.names(short.matrix) == ID(x)))
  {
    x@ps_bg_size <- as.integer(short.matrix[ID(x),"BG_SIZE"])
    x@ps_bg_avg <-  as.numeric(short.matrix[ID(x),"BG_MEAN"])
    x@ps_bg_std_dev <-  as.numeric(short.matrix[ID(x),"BG_STDEV"])
  }
  else
  {
    warning(paste("No background values found for", ID(x), name(x)))
  }

  
  return(x)
})

#' Normalize a `PSMatrix` object
#' 
#' This method normalizes the matrix within a `PSMatrix` object using 
#' column wise normalization and log-transformation. 
#' 
#' @param x An object of `PSMatrix` class.
#' 
#' @return The updated `PSMatrix` object with normalized columns. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' .ps_norm_matrix(pfm1)
#' 
#' @export
#' @importMethodsFrom TFBSTools Matrix
setMethod(".ps_norm_matrix", "PSMatrix", function(x){
  
  mx <- Matrix(x)
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
  mx <- mx + x@.PS_PSEUDOCOUNT
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
  mx <- log(mx)
  
  Matrix(x) <- mx
  
  return(x)
  
})

#' Perform a Scan of DNA Sequences Using a `PSMatrix` Object
#' 
#' This method performs a scan of DNA sequences using a `PSMatrix` object 
#' to identify potential hits based on the matrix score. 
#' 
#' @param x An object of `PSMatrix` class.
#' @param seqs `A DNAstringSet` object representing the sequences to be
#'    scanned.
#' @param BG A logical value indincating whether the scan is for background
#'    sequences. Default is set to `FALSE`.
#' 
#' @return The updated `PSMatrix` object with the scanning results.
#' 
#' @export
setMethod("ps_scan", "PSMatrix", function(x, seqs, BG = FALSE){
  
  if(!is(seqs, "DNAStringSet"))
    stop("seqs is not an object of DNAStringSet class")
  
  rc_x <- reverseComplement(x)
  
  Margs <- list(numx = as.numeric(Matrix(x)), 
                numx_rc = as.numeric(Matrix(rc_x)),
                ncolx = (0:(ncol(Matrix(x)) - 1))*length(.PS_ALPHABET(x)), 
                AB = .PS_ALPHABET(x)) 
  
  if(BG == FALSE)
    x@ps_seq_names <- names(seqs)
  
  seqs <- as.character(seqs)
  
  res <- mapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs)
  
  x <- .ps_add_hits(x, Score = as.numeric(res["score",]), 
                 Strand = as.character(res["strand",]), 
                 Pos = as.integer(res["pos",]), 
                 Oligo = as.character(res["oligo",]), BG = BG)
  
  return(x)
  
})

#' Scan a Single DNA Sequence Using a PSMatrix Object
#' 
#' @param x A `PSMatrix` object containing the position-specific scoring matrix.
#' @param Seq Character. A single DNA sequence represented as a string.
#' @param numx Numeric. Scoring value for the forward strand.
#' @param numx_rc Numeric. Scoring value for the reverse strand. 
#' @param ncolx Integer. Indicates the number of columns 
#' @param AB A named integer vector where names correspond to DNA bases 
#'    and values indicate their respective positions in the scoring matrix.
#' 
#' @return A list containing:
#' \itemize{
#'   \item `score`: the highest score obtained.
#'   \item `strand`: the strand with the highest score.
#'   \item `pos`: the position of the best match in the sequence. 
#'   \item `oligo`: the subsequence corresponding to the best match.}
#' 
#' @importMethodsFrom Biostrings maxScore minScore
#' @export
setMethod(".ps_scan_s", "PSMatrix", function(x, Seq, numx, numx_rc, ncolx, AB){
  
  subS <- strsplit(substring(Seq, seq_len((nchar(Seq) - length(x) + 1)), 
                             length(x):nchar(Seq)),"",
                   fixed = TRUE)
  prot <- numeric(1)
  
  scores <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, 
                   x = numx, AB = AB, ncolx = ncolx)
  scores_rc <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, 
                      x = numx_rc, AB = AB, ncolx = ncolx)
  
  mscore_pos <- which.max(scores)
  mscore_rc_pos <- which.max(scores_rc)
  
  res <- list(score = numeric(), strand = character(), pos = integer(), 
              oligo = character())
  
  if(scores[mscore_pos] >= scores_rc[mscore_rc_pos])
  {
    res$score <- scores[mscore_pos]
    res$strand <- "+"
    res$pos <- mscore_pos
    res$oligo <- paste(subS[[mscore_pos]], collapse = '')
  }
  else
  {
    res$score <- scores_rc[mscore_rc_pos]
    res$strand <- "-"
    res$pos <- mscore_rc_pos
    res$oligo <- paste(subS[[mscore_rc_pos]], collapse = '')
  }
  
  return(res)
})

#' Assign a Score to a DNA Sequence 
#' 
#' This function computes the score for a given DNA subsequence based on 
#' a position-specific scoring matrix.
#' 
#' @param S Character. Represents a DNA sequence.
#' @param x A `PSMatrix` object containing the position-specific scoring matrix.
#' @param AB A named integer vector where names correspond to DNA bases 
#'    and values indicate their respective positions in the scoring matrix.
#' @param ncolx description
#' 
#' @return A single numeric value representing the score of the input 
#'    DNA subsequence.
#' 
#' @importMethodsFrom Biostrings reverseComplement 
#' @export

#setMethod(".ps_assign_score", "PSMatrix", function(x, S){
#  sum(Matrix(x)[matrix(data = c(.PS_ALPHABET(x)[S], 1:length(x)), ncol = 2, nrow = length(x))])
#})

.ps_assign_score <- function(S, x, AB,ncolx)
{
  sum(x[ncolx+AB[S]]) #Assign score to oligo
}

#' Validate a `PSMatrix` object
#' 
#' This function checks if the input `PSMatrix` object properties meet the ones
#' desired.
#' 
#' @param object An object of class `PSMatrix`.
#' 
#' @details
#' The function `validPSMatrix` validates many properties of a `PSMatrix` 
#' object, including: 
#' \itemize{
#'   \item `ps_bg_avg`,`ps_fg_avg`, and `ps_bg_std_dev` must be of length 1.
#'   \item `ps_bg_avg` and `ps_bg_std_dev` value must be between 0 and 1, 
#'   excluding 0 for `ps_bg_std_dev`.
#'   \item the length of `ps_hits_pos`, `ps_hits_strand`, and `ps_hits_score` 
#'   vectors must be equal.}
#' 
#' @return If all the checks are satisfied, returns `TRUE`. Otherwise, a
#' string describing the reason of failure. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' validPSMatrix(pfm1)
#' 
#' @export
validPSMatrix <- function(object)
{
  if(length(object@ps_bg_avg) != 1)
    return("Background average must be of length 1")
  if(length(object@ps_fg_avg) != 1)
    return("Foreground average must be of length 1")
  if(length(object@ps_bg_std_dev) != 1)
    return("Background stdev must be of length 1")
  if((object@ps_bg_avg < 0 || object@ps_bg_avg > 1) && !is.na(object@ps_bg_avg)) 
    return(paste("Invalid value for Background average: ", object@ps_bg_avg))
  if((object@ps_bg_std_dev <= 0 || object@ps_bg_std_dev > 1) && !is.na(object@ps_bg_std_dev))
    return(paste("Invalid value for Background stddev: ", object@ps_bg_std_dev))
 # if(object@ps_bg_size < 1000 && !is.na(object@ps_bg_size))
  #  return(paste("Invalid value for Background size: ", object@ps_bg_size, " Background must be of at least 1000 sequences"))
  if(length(object@ps_hits_pos) != length(object@ps_hits_strand) || length(object@ps_hits_pos) != length(object@ps_hits_score))
    return(paste("Invalid PSMatrix object: different values for hits, strands and scores vectors"))
  
  TRUE
}

#' @export
setValidity("PSMatrix", validPSMatrix)

#' Display Details of a `PSMatrix` object
#' 
#' This method displays information about a `PSMatrix` object, including 
#' background and foreground statistics, satndard deviation, background 
#' and foreground size, z-scores and P-Value. 
#' 
#' @param object An object of class `PSMatrix`.
#' 
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_fg_avg}}, 
#'    \code{\link{ps_bg_std_dev}}, \code{\link{ps_bg_size}}, 
#'    \code{\link{ps_fg_size}}, \code{\link{ps_zscore}}, \code{\link{ps_pvalue}}
#'    
#' @return An object of class `PSMatrix` with its details displyed. 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' show(pfm1)
#' 
#' @export
#' @importMethodsFrom methods show
setMethod("show", "PSMatrix", function(object) {
  
  callNextMethod()
  
  cat(
      "\nPscan Background Average: ", ps_bg_avg(object),
      "\nPscan Foreground (your sample) Average: ", ps_fg_avg(object),
      "\nPscan Backgroun Stdev: ", ps_bg_std_dev(object),
      "\nPscan Background Size: ", ps_bg_size(object),
      "\nPscan Foreground (your sample) Size: ", ps_fg_size(object),
      "\nPscan Zscore: ", ps_zscore(object),
      "\nPscan pvalue (z-test): ", ps_pvalue(object),
      sep = ""
  )
})

#' Generic Setter Method for Background Mean Value of an object
#' 
#' @param x An object for which the background average is to be set.
#' @param ... Additional arguments.
#' @param value The value for the background average. 
#' 
#' @return An object of the same class of the input with the modified 
#'    background average value. 
#'    
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_avg<-`(pfm1, value = 0.73473) 
#' 
#' @export
setGeneric("ps_bg_avg<-", function(x, ..., value) standardGeneric("ps_bg_avg<-"))

#' Generic Setter Method for Background Standard Deviation Value of an object
#' 
#' @param x An object for which the background standard deviation
#'    is to be set.
#' @param ... Additional arguments.
#' @param value The value for the background standard deviation. 
#' 
#' @return An object of the same class of the input with the modified 
#'    background standard deviation value. 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_std_dev<-`(pfm1, value = 0.08568925)
#' 
#' @export
setGeneric("ps_bg_std_dev<-", function(x, ..., value) standardGeneric("ps_bg_std_dev<-"))

#' Generic Setter Method for Background Size Value of an object
#' 
#' @param x An object for which the background size value is to be set.
#' @param ... Additional arguments.
#' @param value The value for the background size. 
#' 
#' @return An object of the same class of the input with the modified 
#'    background size. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_size<-`(pfm1, value = 25629)
#' 
#' 
#' @export
setGeneric("ps_bg_size<-", function(x, ..., value) standardGeneric("ps_bg_size<-"))

#' Setter Method for Background Average in PSMatrix
#' 
#' This method specifically set the background average for an object of 
#' class `PSMatrix`.
#' 
#' @param x An object of class `PSMatrix`.
#' @param value Numeric. The mean background value to be set.
#' 
#' @return The modified `PSMatrix` object with the updated background average. 
#' 
#' @seealso \code{\link{ps_bg_avg}} to retrieve the background average value,
#'    and \code{\link{validObject}} to see which checks are performed on the
#'    input `PSMatrix` object.
#'    
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_avg<-`(pfm1, value = 0.73473) 
#' 
#' @export
setReplaceMethod("ps_bg_avg", "PSMatrix", function(x,value){
  
  x@ps_bg_avg <- value
  validObject(x)
  x
})

#' Setter Method for Background Standard Deviation in PSMatrix
#' 
#' This method specifically set the background standard deviation
#' for an object of class `PSMatrix`.
#' 
#' @param x An object of class `PSMatrix`.
#' @param value Numeric. The standard deviation background value to be set.
#' 
#' @return The modified `PSMatrix` object with the updated background 
#'    standard deviation. 
#' 
#' @seealso \code{\link{ps_bg_std_dev}} to retrieve the background 
#'    standard deviation, and \code{\link{validObject}} to see which checks 
#'    are performed on the input `PSMatrix` object.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_std_dev<-`(pfm1, value = 0.08568925)
#' 
#' @export
setReplaceMethod("ps_bg_std_dev", "PSMatrix", function(x,value){
  
  x@ps_bg_std_dev <- value
  validObject(x)
  x
})

#' Setter Method for Background Size in PSMatrix
#' 
#' This method specifically set the background size for an object 
#' of class `PSMatrix`.
#' 
#' @param x An object of class `PSMatrix`.
#' @param value Numeric. The size background value to be set.
#' 
#' @return The modified `PSMatrix` object with the updated background size. 
#' 
#' @seealso \code{\link{ps_bg_size}} to retrieve the background size, 
#'    and \code{\link{validObject}} to see which checks are performed on the
#'    input `PSMatrix` object.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rda", package = "PscanR")
#' load(pfm1_path)
#' `ps_bg_size<-`(pfm1, value = 25629)
#' 
#' @export
setReplaceMethod("ps_bg_size", "PSMatrix", function(x,value){
  
  x@ps_bg_size <- value
  validObject(x)
  x
})

#' Coercion from PFMatrix to PSMatrix
#' 
#' @param from A `PFMatrix` object to be converted to a `PSMatrix` one.
#' 
#' @return A `PSMatrix` object created from the input `PFMatrix` object.
#' 
#' @export
#' @importFrom TFBSTools PFMatrix
setAs("PFMatrix", "PSMatrix", function(from){
  
  #.ps_norm_matrix(new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
  #                    ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01))
  
  PSMatrix(from)
  
})

#PSMatrix <- function(pfm, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
#                     ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01, ...)

#' Create a PSMatrixList object
#' 
#' This method creates a `PSMatrixList` object starting from different 
#' `PFMatrix` object.
#' 
#' @param from An object of class `PFMatrix`.
#' 
#' @return An object of class `PSMatrixList`.
#' 
#' @export
#' @importFrom TFBSTools PFMatrixList
setAs("PFMatrixList", "PSMatrixList", function(from){
  
  to <- lapply(from, as, "PSMatrix")
  
  do.call(PSMatrixList, to)
  
})

