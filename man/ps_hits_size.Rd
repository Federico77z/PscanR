% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{ps_hits_size}
\alias{ps_hits_size}
\title{Retrieve the Hits Size Value for Motif Enrichment Analysis 
(Generic Function)}
\usage{
ps_hits_size(x, ...)
}
\arguments{
\item{x}{An object containing motif scan results.}

\item{...}{Additional arguments passed to specific methods.}
}
\value{
\itemize{
  \item If `x` is a `PSMatrix` object (Pscan result), an **integer value**.
  \item If `x` is another supported class, the return format may differ.
}
}
\description{
`ps_hits_size` is a **generic function** that retrieves the number of hits
for motif enrichment analysis. Methods should be 
implemented for specific object classes that store scan results.
}
\details{
The hits size value represent the total number of motif hits detected in 
the input promoter sequences.
}
\examples{
pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
pfm1 <- readRDS(pfm1_path)
ps_hits_size(pfm1)

}
