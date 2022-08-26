#' @export
#' @import methods
#' @importClassesFrom TFBSTools TFBSTools

.PSMatrix <- setClass("PSMatrix", slots = representation(ps_bg_avg="numeric", 
                                                         ps_bg_std_err="numeric", 
                                                         ps_bg_size="integer", 
                                                         ps_hits_pos="integer", 
                                                         ps_hits_strand="character", 
                                                         ps_hits_score="numeric",
                                                         .PS_PSEUDOCOUNT="numeric",
                                                         .PS_ALPHABET="integer"), 
                      contains="PFMatrix")

PSMatrix <- function(ps_bg_avg = NA, ps_bg_std_err = NA, ps_bg_size = NA, 
                     .PS_PSEUDOCOUNT = 0.01, .PS_ALPHABET = setNames(1:4, c("A","C","G","T")), ...)
{
  pfm <- PFMatrix(...)
  .PSMatrix(pfm, ps_bg_avg = ps_bg_avg, ps_bg_std_err = ps_bg_std_err, ps_bg_size = ps_bg_size, 
            .PS_PSEUDOCOUNT = .PS_PSEUDOCOUNT, .PS_ALPHABET=.PS_ALPHABET)
}

#' @export

setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' @export
setGeneric("ps_bg_std_err", function(x, ...) standardGeneric("ps_bg_std_err"))

#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

#' @export
setGeneric("ps_hits_size", function(x, ...) standardGeneric("ps_hits_size"))

#' @export
setGeneric(".PS_PSEUDOCOUNT", function(x, ...) standardGeneric(".PS_PSEUDOCOUNT"))

#' @export
setGeneric(".PS_ALPHABET", function(x, ...) standardGeneric(".PS_ALPHABET"))

#' @export
setGeneric(".ps_norm_matrix", function(x, ...) standardGeneric(".ps_norm_matrix"))

#' @export
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))

#' @export
setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))


#setGeneric(".ps_assign_score", function(x, ...) standardGeneric(".ps_assign_score"))

#' @export
setGeneric(".ps_add_hit", function(x, ...) standardGeneric(".ps_add_hit"))

#' @export
setMethod("ps_bg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_avg
  
  return(out)
})

#' @export
#' 
setMethod("ps_bg_std_err", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_std_err
  
  return(out)
})

#' @export
#' 
setMethod("ps_bg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_size
  
  return(out)
})

#' @export

setMethod("ps_hits_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- length(x@ps_hits_pos)
  
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
setMethod(".ps_add_hit", "PSMatrix", function(x, Pos, Strand, Score, withDimnames = TRUE) {
  
  x@ps_hits_pos <- Pos
  x@ps_hits_strand <- Strand
  x@ps_hits_score <- Score
  
  return(x)
})

#' @export
#' @importMethodsFrom TFBSTools Matrix

setMethod(".ps_norm_matrix", "PSMatrix", function(x){
  
  mx <- TFBSTools::Matrix(x)
  
  sums <- apply(mx, 2, sum)
  
  mx <- mx / sums
  
  mx <- mx + x@.PS_PSEUDOCOUNT
  
  sums <- apply(mx, 2, sum)
  
  mx <- log(mx / sums)
  
  Matrix(x) <- mx
  
  return(x)
  
})

#' @export

setMethod("ps_scan", "PSMatrix", function(x, seqs){
  
  if(!is(seqs, "DNAStringSet"))
    stop("seqs is not an object of DNAStringSet class")
  
  rc_x <- reverseComplement(x)
  
  Margs = list(numx = as.numeric(Matrix(x)), numx_rc = as.numeric(Matrix(rc_x)),
               ncolx = (0:(ncol(Matrix(x)) - 1))*length(.PS_ALPHABET(x)), AB = .PS_ALPHABET(x)) 
  
  seqs <- as.character(seqs)
  
  res <- mapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs)
  #bpmapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs, BPPARAM = MulticoreParam())
  
  x <- .ps_add_hit(x, Score = as.numeric(res["score",]), 
                 Strand = as.character(res["strand",]), 
                 Pos = as.integer(res["pos",]))
  
  #.ps_add_hit(x, Score = 5, 
   #           Strand = "ciao",
  #           Pos = 3L)

 # print("qui")
  
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
  sum(x[ncolx+AB[S]]) #fastest way to assign score to an oligo I was able to figure out (without recurring to C implementation)
}

#' @export

validPSMatrix <- function(object)
{
  if(length(object@ps_bg_avg) != 1)
    return("Background average must be of length 1")
  if(length(object@ps_bg_std_err) != 1)
    return("Background stderr must be of length 1")
  if(length(object@ps_bg_std_err) != 1)
    return("Background size must be of length 1")
  if((object@ps_bg_avg < 0 | object@ps_bg_avg > 1) & !is.na(object@ps_bg_avg)) 
    return(paste("Invalid value for Background average: ", object@ps_bg_avg))
  if((object@ps_bg_std_err < 0 | object@ps_bg_std_err > 1) & !is.na(object@ps_bg_std_err))
    return(paste("Invalid value for Background stderr: ", object@ps_bg_std_err))
  if(object@ps_bg_size <= 1000 & !is.na(object@ps_bg_size))
    return(paste("Invalid value for Background size: ", object@ps_bg_size, " Background must be of at least 1000 sequences"))
  if(length(object@ps_hits_pos) != length(object@ps_hits_strand) | length(object@ps_hits_pos) != length(object@ps_hits_score))
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
      "\nPscan Background Average: ", ps_bg_avg(object), "\n",
      "\nPscan Backgroun StdErr: ", ps_bg_std_err(object), "\n",
      "\nPscan Background Size: ", ps_bg_size(object), "\n",
      "\nPscan Hits Size: ", ps_hits_size(object), "\n",
      sep = ""
  )
})

#' @export
setGeneric(".ps_bg_avg<-", function(x, ..., value) standardGeneric(".ps_bg_avg<-"))

#' @export
setGeneric(".ps_bg_std_err<-", function(x, ..., value) standardGeneric(".ps_bg_std_err<-"))

#' @export
setGeneric(".ps_bg_size<-", function(x, ..., value) standardGeneric(".ps_bg_size<-"))

#' @export

setReplaceMethod(".ps_bg_avg", "PSMatrix", function(x,value){
  
  x@ps_bg_avg <- value
  validObject(x)
  x
})

#' @export

setReplaceMethod(".ps_bg_std_err", "PSMatrix", function(x,value){
  
  x@ps_bg_std_err <- value
  validObject(x)
  x
})

#' @export

setReplaceMethod(".ps_bg_size", "PSMatrix", function(x,value){
  
  x@ps_bg_size <- value
  validObject(x)
  x
})

#' @exportMethods coerce

setAs("PFMatrix", "PSMatrix", function(from){
  
  .ps_norm_matrix(new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_bg_std_err = as.numeric(NA), 
                      ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = as.numeric(0.01), 
                      .PS_ALPHABET = setNames(1:4, c("A","C","G","T"))))
  
})

