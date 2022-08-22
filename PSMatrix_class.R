#' @export
#' @import methods
#' @importClassesFrom TFBSTools TFBSTools

.PSMatrix <- setClass("PSMatrix", slots = representation(ps_bg_avg="numeric", 
                                                         ps_bg_std_err="numeric", 
                                                         ps_bg_size="integer", 
                                                         ps_hits_pos="integer", 
                                                         ps_hits_strand="character", 
                                                         ps_hits_score="numeric",
                                                         .PS_PSEUDOCOUNT="numeric"), 
                      contains="PFMatrix")

PSMatrix <- function(ps_bg_avg = NA, ps_bg_std_err = NA, ps_bg_size = NA, .PS_PSEUDOCOUNT = 0.01, ...)
{
  pfm <- PFMatrix(...)
  .PSMatrix(pfm, ps_bg_avg = ps_bg_avg, ps_bg_std_err = ps_bg_std_err, ps_bg_size = ps_bg_size, .PS_PSEUDOCOUNT = .PS_PSEUDOCOUNT)
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
setGeneric(".ps_norm_matrix", function(x, ...) standardGeneric(".ps_norm_matrix"))

#' @export
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))

#' @export
setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))

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
  
  return(x)
})

#' @importMethodsFrom Biostrings maxScore minScore
#' @export

setMethod(".ps_scan_s", "PSMatrix", function(x, Seq){
  
  from <- 1
  to <- length(Seq) - length(x)
  width <- length(x)
  
  vapply(from:to, )
  
  return(x)
})

#' @importMethodsFrom Biostrings reverseComplement 
#' @export

setMethod(".ps_assign_score", "PSMatrix", function(x, Seq){
  
  sum()
  
})

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
                      ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = as.numeric(0.01)))
  
})

