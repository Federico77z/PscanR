#' @export
#' @import methods
#' @importClassesFrom TFBSTools TFBSTools

.PSMatrix <- setClass("PSMatrix", slots = representation(ps_bg_avg="numeric", ps_bg_std_err="numeric", ps_bg_size="integer"), contains="PFMatrix")

PSMatrix <- function(ps_bg_avg = NA, ps_bg_std_err = NA, ps_bg_size = NA, ...)
{
  pfm <- PFMatrix(...)
  .PSMatrix(pwm, ps_bg_avg = ps_bg_avg, ps_bg_std_err = ps_bg_std_err, ps_bg_size = ps_bg_size)
}

#' @export

setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' @export
setGeneric("ps_bg_std_err", function(x, ...) standardGeneric("ps_bg_std_err"))

#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

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
  
  new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_bg_std_err = as.numeric(NA), ps_bg_size = as.integer(NA))
  
})