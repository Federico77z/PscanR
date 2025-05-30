% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-PSMatrix.R
\name{ps_fg_size,PSMatrix-method}
\alias{ps_fg_size,PSMatrix-method}
\title{Get Foreground Size}
\usage{
\S4method{ps_fg_size}{PSMatrix}(x, withDimnames = TRUE)
}
\arguments{
\item{x}{A `PSMatrix` object.}

\item{withDimnames}{Logical, whether to include dimension names in the 
output, if they exist in the object.
Default set to `TRUE`.}
}
\value{
An integer representing the total number of promoter sequences used 
   as input to Pscan.
}
\description{
Retrieves the number of promoter sequences given as input to the Pscan 
algorithm in a `PSMatrix` object. This represents the total number of 
sequences analyzed for motif enrichment.
}
\examples{
pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
pfm1 <- readRDS(pfm1_path)
ps_fg_size(pfm1)

}
