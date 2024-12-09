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
#' @slot ps_hits_scroe Numeric. Score information of hits. 
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
#' @export
setGeneric("ps_zscore", function(x, ...) standardGeneric("ps_zscore"))

#' @describeIn ps_generics 
#' Retrieve the P-Value from an object
#'
#' @return For `ps_pvalue`: a numeric value representing the P-Value.
#' @export
setGeneric("ps_pvalue", function(x, ...) standardGeneric("ps_pvalue"))

#' @describeIn ps_generics 
#' Retrieve the average background value from an object 
#'
#' @return For `ps_bg_avg`: a numeric value representing the average background.
#' @export
setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' @describeIn ps_generics 
#' Retrieve the average foreground value from an object
#'
#' @return For `ps_fg_avg`: a numeric value representing average foreground.
#' @export
setGeneric("ps_fg_avg", function(x, ...) standardGeneric("ps_fg_avg"))

#' @describeIn ps_generics 
#' Retrieve the background standard deviation from an object
#'
#' @return For `ps_bg_std_dev`: a numeric value representing the background 
#'    standard deviation.
#' @export
setGeneric("ps_bg_std_dev", function(x, ...) standardGeneric("ps_bg_std_dev"))

#' @describeIn ps_generics 
#' Retrieve the background size from an object
#'
#' @return For `ps_bg_size`: an integer value representing the size of the 
#'    background dataset.
#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

#' @describeIn ps_generics 
#' Retrieve the foreground size from an object
#'
#' @return For `ps_fg_size`: an integer value representing the size pf 
#'    the foreground dataset.
#' @export
setGeneric("ps_fg_size", function(x, ...) standardGeneric("ps_fg_size"))

#' @describeIn ps_generics 
#' Retrieve the hits size from an object
#'
#' @return For `ps_hits_size`: an integer value representing the hits size.
#' @export
setGeneric("ps_hits_size", function(x, ...) standardGeneric("ps_hits_size"))

#' @describeIn ps_generics 
#' Retrieve the hits score from an object
#'
#' @return For `ps_hits_score`: a numeric vector representing the score for 
#'    each hit.
#' @export
setGeneric("ps_hits_score", function(x, ...) standardGeneric("ps_hits_score"))

#' @describeIn ps_generics 
#'
#' @return For `ps_hits_z`: a numeric vector.
#' @export
setGeneric("ps_hits_z", function(x, ...) standardGeneric("ps_hits_z"))

#' @describeIn ps_generics 
#' Retrieve the strand of hits from an object
#'
#' @return For `ps_hits_strand`: a character vector representing the strand 
#'    of each hit.
#' @export
setGeneric("ps_hits_strand", function(x, ...) standardGeneric("ps_hits_strand"))

#' @describeIn ps_generics 
#' Retrieve the position of hits from an object
#'
#' @return For `ps_hits_pos`: an integer vector representing the position of
#'    each hit.
#' @export
setGeneric("ps_hits_pos", function(x, ...) standardGeneric("ps_hits_pos"))

#' @describeIn ps_generics 
#' Retrieve oligonucleotide sequence of hits from an object
#'
#' @return For `ps_hits_oligo`: a character vector representing the 
#'    oligonucleotide sequence of each hit. 
#' @export
setGeneric("ps_hits_oligo", function(x, ...) standardGeneric("ps_hits_oligo"))

#' @describeIn ps_generics 
#' Retrieve the hits table from an object
#'
#' @return For `ps_hits_table`: a `data.frame` of hits.
#' @export
setGeneric("ps_hits_table", function(x, ...) standardGeneric("ps_hits_table"))
# not sure

#' @describeIn ps_generics 
#' Retrieve the sequence name from an object
#'
#' @return For `ps_seq_names`: a character vector of sequence names.
#' @export
setGeneric("ps_seq_names", function(x, ...) standardGeneric("ps_seq_names"))

#' @describeIn ps_generics
#' Retrieve the pseudocount value
#' 
#' @return For `.PS_PSEUDOCOUNT`: a numeric value representing the pseudocount.
#' @export
setGeneric(".PS_PSEUDOCOUNT", function(x, ...) standardGeneric(".PS_PSEUDOCOUNT"))

#' @describeIn ps_generics
#' Retrieve the numbers representing the DNA alphabet
#'  
#' @return For `.PS_ALPHABET`: integer values.
#' @export
setGeneric(".PS_ALPHABET", function(x, ...) standardGeneric(".PS_ALPHABET"))

#' Normalize an Object 
#' 
#' @param x `PSMatrix` Object 
#' @param ... Additional parameters
#' 
#' @export
setGeneric(".ps_norm_matrix", function(x, ...) standardGeneric(".ps_norm_matrix"))

#' @export
setGeneric(".ps_seq_names", function(x, out) standardGeneric(".ps_seq_names"))

#' @describeIn ps_generics 
#' Perform a scan operation on an Object 
#'
#' @return A `data.frame` of hits. 
#' @export
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))

#' Retrieve the background values from a table
#' 
#' @param x A `data.frame`
#' @param ... Additional parameters
#' 
#' @export
setGeneric(".ps_bg_from_table", function(x, ...) standardGeneric(".ps_bg_from_table"))

#' @export
setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))

#' @export
setGeneric(".ps_norm_score", function(x, ...) standardGeneric(".ps_norm_score"))

#setGeneric(".ps_assign_score", function(x, ...) standardGeneric(".ps_assign_score"))

#' @export
setGeneric(".ps_add_hits", function(x, ...) standardGeneric(".ps_add_hits"))

#' Get Background Average Score
#'
#' Retrieves the background average score stored in a `PSMatrix` object.
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output,
#'    if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the background average score. 
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
#' @param withDimnames Logical, whether to include dimension names in the output,
#'    if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the foreground average score. 
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
#' @param withDimnames Logical, whether to include dimension names in the output,
#'    if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the `PSMatrix` Z-score.
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
#' @param withDimnames Logical, whether to include dimension names in the output,
#'    if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the `PSMatrix` P-Value.
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
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
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
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric value representing the `PSMatrix` background standard deviation.
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
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric value representing the `PSMatrix` background size.
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
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric value representing the `PSMatrix` foreground size.
#'
#' @export
setMethod("ps_fg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_size
  
  return(out)
})

#' @export
setMethod("ps_hits_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- length(x@ps_hits_pos)
  
  return(out)
})

#' @export
setMethod("ps_hits_score", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_score
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' @export
setMethod("ps_hits_z", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- (x@ps_hits_score - x@ps_bg_avg) / x@ps_bg_std_dev
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' @export
setMethod("ps_hits_strand", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_strand
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' @export
setMethod("ps_hits_pos", "PSMatrix", function(x, pos_shift = 0L, 
                                              withDimnames = TRUE) {
  out <- x@ps_hits_pos + pos_shift
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' @export
setMethod("ps_seq_names", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_seq_names
  
  return(out)
})

#' @export

setMethod(".PS_PSEUDOCOUNT", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_PSEUDOCOUNT
  
  return(out)
})

#' @export

setMethod(".PS_ALPHABET", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_ALPHABET
  
  return(out)
})

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

#' @export
setMethod(".ps_norm_score", "PSMatrix", function(x) {
  
  ps_score <- 1 + ((maxScore(Matrix(x)) - ps_hits_score(x)) / 
         (minScore(Matrix(x)) - maxScore(Matrix(x))))
  
  return(ps_score)
})

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

#' @importMethodsFrom Biostrings reverseComplement 
#' @export

#setMethod(".ps_assign_score", "PSMatrix", function(x, S){
#  sum(Matrix(x)[matrix(data = c(.PS_ALPHABET(x)[S], 1:length(x)), ncol = 2, nrow = length(x))])
#})

.ps_assign_score <- function(S, x, AB,ncolx)
{
  sum(x[ncolx+AB[S]]) #Assign score to oligo
}

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

#' @export

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

#' @export
setGeneric("ps_bg_avg<-", function(x, ..., value) standardGeneric("ps_bg_avg<-"))

#' @export
setGeneric("ps_bg_std_dev<-", function(x, ..., value) standardGeneric("ps_bg_std_dev<-"))

#' @export
setGeneric("ps_bg_size<-", function(x, ..., value) standardGeneric("ps_bg_size<-"))

#' @export

setReplaceMethod("ps_bg_avg", "PSMatrix", function(x,value){
  
  x@ps_bg_avg <- value
  validObject(x)
  x
})

#' @export

setReplaceMethod("ps_bg_std_dev", "PSMatrix", function(x,value){
  
  x@ps_bg_std_dev <- value
  validObject(x)
  x
})

#' @export

setReplaceMethod("ps_bg_size", "PSMatrix", function(x,value){
  
  x@ps_bg_size <- value
  validObject(x)
  x
})

#' @export

setAs("PFMatrix", "PSMatrix", function(from){
  
  #.ps_norm_matrix(new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
  #                    ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01))
  
  PSMatrix(from)
  
})

#PSMatrix <- function(pfm, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
#                     ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01, ...)

#' @export

setAs("PFMatrixList", "PSMatrixList", function(from){
  
  to <- lapply(from, as, "PSMatrix")
  
  do.call(PSMatrixList, to)
  
})

