#' @export
#' @import methods
#' @importClassesFrom TFBSTools TFBSTools

.PSMatrix <- setClass("PSMatrix", slots = representation(ps_bg_avg="numeric",
                                                         ps_fg_avg="numeric",
                                                         ps_bg_std_dev="numeric", 
                                                         ps_bg_size="integer",
                                                         ps_fg_size="integer",
                                                         ps_hits_pos="integer", 
                                                         ps_hits_strand="character", 
                                                         ps_hits_score="numeric",
                                                         ps_zscore="numeric",
                                                         ps_pvalue="numeric",
                                                         ps_seq_names="character",
                                                         .PS_PSEUDOCOUNT="numeric",
                                                         .PS_ALPHABET="integer"), 
                      contains="PFMatrix")

PSMatrix <- function(ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
                     ps_bg_size = as.integer(NA), ps_fg_size = as.integer(NA), ps_zscore = as.numeric(NA),
                     ps_pvalue = as.numeric(NA), ps_seq_names = as.character(NA),
                     .PS_PSEUDOCOUNT = 0.01, .PS_ALPHABET = setNames(1:4, c("A","C","G","T")), ...)
{
  pfm <- PFMatrix(...)
  .PSMatrix(pfm, ps_bg_avg = ps_bg_avg, ps_fg_avg = ps_fg_avg, ps_bg_std_dev = ps_bg_std_dev, ps_bg_size = ps_bg_size, 
            ps_fg_size = ps_fg_size, ps_zscore = ps_zscore, ps_pvalue = ps_pvalue, ps_seq_names = ps_seq_names,
            .PS_PSEUDOCOUNT = .PS_PSEUDOCOUNT, .PS_ALPHABET=.PS_ALPHABET, ps_hits_pos = integer(), 
            ps_hits_strand = character(), ps_hits_score = numeric())
}

.PSMatrixList <-setClass("PSMatrixList", contains ="PFMatrixList")

PSMatrixList <- function(..., use.names = TRUE)
{
  listData = list(...)
  XMatrixList(listData, use.names = use.names, type = "PSMatrixList", matrixClass = "PSMatrix")
}

#' @export
setGeneric("ps_zscore", function(x, ...) standardGeneric("ps_zscore"))

#' @export
setGeneric("ps_pvalue", function(x, ...) standardGeneric("ps_pvalue"))

#' @export
setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' @export
setGeneric("ps_fg_avg", function(x, ...) standardGeneric("ps_fg_avg"))

#' @export
setGeneric("ps_bg_std_dev", function(x, ...) standardGeneric("ps_bg_std_dev"))

#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

#' @export
setGeneric("ps_fg_size", function(x, ...) standardGeneric("ps_fg_size"))

#' @export
setGeneric("ps_hits_size", function(x, ...) standardGeneric("ps_hits_size"))

#' @export
setGeneric("ps_hits_score", function(x, ...) standardGeneric("ps_hits_score"))

#' @export
setGeneric("ps_hits_strand", function(x, ...) standardGeneric("ps_hits_strand"))

#' @export
setGeneric("ps_hits_pos", function(x, ...) standardGeneric("ps_hits_pos"))

#' @export
setGeneric("ps_hits_table", function(x, ...) standardGeneric("ps_hits_table"))

#' @export
setGeneric("ps_seq_names", function(x, ...) standardGeneric("ps_seq_names"))

#' @export
setGeneric(".PS_PSEUDOCOUNT", function(x, ...) standardGeneric(".PS_PSEUDOCOUNT"))

#' @export
setGeneric(".PS_ALPHABET", function(x, ...) standardGeneric(".PS_ALPHABET"))

#' @export
setGeneric(".ps_norm_matrix", function(x, ...) standardGeneric(".ps_norm_matrix"))

#' @export
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))

#' @export
setGeneric(".ps_bg_from_table", function(x, ...) standardGeneric(".ps_bg_from_table"))

#' @export
setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))

#' @export
setGeneric(".ps_norm_score", function(x, ...) standardGeneric(".ps_norm_score"))

#setGeneric(".ps_assign_score", function(x, ...) standardGeneric(".ps_assign_score"))

#' @export
setGeneric(".ps_add_hits", function(x, ...) standardGeneric(".ps_add_hits"))

#' @export
setMethod("ps_bg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_avg
  
  return(out)
})

#' @export
setMethod("ps_fg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_avg
  
  return(out)
})

#' @export
setMethod("ps_zscore", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_zscore
  
  return(out)
})

#' @export
setMethod("ps_pvalue", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_pvalue
  
  return(out)
})

#' @export
#' 
setMethod("ps_bg_std_dev", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_std_dev
  
  return(out)
})

#' @export
#' 
setMethod("ps_bg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_size
  
  return(out)
})

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
  
  return(out)
})

#' @export
setMethod("ps_hits_strand", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_strand
  
  return(out)
})

#' @export
setMethod("ps_hits_pos", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_pos
  
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

setMethod("ps_hits_table", "PSMatrix", function(x, pos_shift = 0L, withDimnames = TRUE) {
  
  out <- data.frame("SCORE" = x@ps_hits_score, "POS" = x@ps_hits_pos + as.integer(pos_shift), "STRAND" = x@ps_hits_strand,
                    row.names = x@ps_seq_names)
  
  out <- out[with(out, order(SCORE, POS, decreasing = c(TRUE,FALSE))),]
  
  return(out)
})

#' @export
setMethod(".ps_add_hits", "PSMatrix", function(x, Pos, Strand, Score, BG = FALSE, withDimnames = TRUE) {
  
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
      ztest <- z.test(x@ps_hits_score, mu = x@ps_bg_avg, sigma.x = x@ps_bg_std_dev, alternative = "greater")
      
      x@ps_zscore <- ztest$statistic["z"]
      x@ps_pvalue <- as.numeric(ztest$p.value)
      x@ps_fg_avg <- mean(x@ps_hits_score, na.rm = TRUE)
      x@ps_fg_size <- length(x@ps_hits_pos)
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
  
  if(any(row.names(short.matrix) == name(x)))
  {
    x@ps_bg_size <- as.integer(short.matrix[name(x),"BG_SIZE"])
    x@ps_bg_avg <-  as.numeric(short.matrix[name(x),"BG_MEAN"])
    x@ps_bg_std_dev <-  as.numeric(short.matrix[name(x),"BG_STDEV"])
  }
  else
  {
    warning(paste("No background values found for", name(x), ID(x)))
  }

  
  return(x)
})

#' @export
#' @importMethodsFrom TFBSTools Matrix

setMethod(".ps_norm_matrix", "PSMatrix", function(x){
  
  mx <- Matrix(x)
  
  #sums <- apply(mx, 2, sum)
  
  #mx <- mx / sums
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
  mx <- mx + x@.PS_PSEUDOCOUNT
  
#  sums <- apply(mx, 2, sum)
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
#  mx <- log(mx / sums)
  
  mx <- log(mx)
  
  Matrix(x) <- mx
  
  return(x)
  
})

#' @export

setMethod("ps_scan", "PSMatrix", function(x, seqs, BG = FALSE){
  
  if(!is(seqs, "DNAStringSet"))
    stop("seqs is not an object of DNAStringSet class")
  
  rc_x <- reverseComplement(x)
  
  Margs = list(numx = as.numeric(Matrix(x)), numx_rc = as.numeric(Matrix(rc_x)),
               ncolx = (0:(ncol(Matrix(x)) - 1))*length(.PS_ALPHABET(x)), AB = .PS_ALPHABET(x)) 
  
  if(BG == FALSE)
    x@ps_seq_names = names(seqs)
  
  seqs <- as.character(seqs)
  
  res <- mapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs)
  
  x <- .ps_add_hits(x, Score = as.numeric(res["score",]), 
                 Strand = as.character(res["strand",]), 
                 Pos = as.integer(res["pos",]), BG = BG)
  
  return(x)
  
})

#' @importMethodsFrom Biostrings maxScore minScore
#' @export

setMethod(".ps_scan_s", "PSMatrix", function(x, Seq, numx, numx_rc, ncolx, AB){
  
  subS <- strsplit(substring(Seq, 1:(nchar(Seq) - length(x) + 1), length(x):nchar(Seq)),"",
                                      fixed = TRUE)
  prot <- numeric(1)
  
  scores <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, x = numx, AB = AB, ncolx = ncolx)
  scores_rc <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, x = numx_rc, AB = AB, ncolx = ncolx)
  
  mscore_pos <- which.max(scores)
  mscore_rc_pos <- which.max(scores_rc)
  
  res <- list(score = numeric(), strand = character(), pos = integer())
  
  if(scores[mscore_pos] >= scores_rc[mscore_rc_pos])
  {
    res$score <- scores[mscore_pos]
    res$strand <- "+"
    res$pos <- mscore_pos
  }
  else
  {
    res$score <- scores_rc[mscore_rc_pos]
    res$strand <- "-"
    res$pos <- mscore_rc_pos
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
  if(object@ps_bg_size < 1000 && !is.na(object@ps_bg_size))
    return(paste("Invalid value for Background size: ", object@ps_bg_size, " Background must be of at least 1000 sequences"))
  if(length(object@ps_hits_pos) != length(object@ps_hits_strand) || length(object@ps_hits_pos) != length(object@ps_hits_score))
    return(paste("Invalid PSMatrix object: different values for hits, strands and scores vectors"))
  
  TRUE
}

#' @export
setValidity("PSMatrix", validPSMatrix)

#' @export
#' @importMethodsFrom PFMatrix show

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

#' @exportMethods coerce

setAs("PFMatrix", "PSMatrix", function(from){
  
  .ps_norm_matrix(new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
                      ps_zscore = as.numeric(NA), ps_pvalue = as.numeric(NA), ps_bg_size = as.integer(NA), 
                      ps_fg_size = as.integer(NA), ps_seq_names = as.character(NA),
                      .PS_PSEUDOCOUNT = as.numeric(0.01), .PS_ALPHABET = setNames(1:4, c("A","C","G","T"))))
  
})

#' @exportMethods coerce

setAs("PFMatrixList", "PSMatrixList", function(from){
  
  to <- lapply(from, as, "PSMatrix")
  
  do.call(PSMatrixList, to)
  
})

