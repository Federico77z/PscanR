% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{ps_hits_pos_bg}
\alias{ps_hits_pos_bg}
\title{Retrieve Motif Hit Positions on Background Dataset (Generic Function)}
\usage{
ps_hits_pos_bg(x, ...)
}
\arguments{
\item{x}{An object containing motif scan results.}

\item{...}{Additional arguments passed to specific methods.}
}
\value{
\itemize{
  \item If `x` is a `PSMatrix` object: an **integer vector**, where 
  names correspond to promoter sequence identifiers, and values represent 
  motif hit positions along the background dataset.
  \item If `x` is another supported class, the return format may differ.
}
}
\description{
`ps_hits_pos_bg` is a **generic function** that retrieves the positions of 
motif hits in a given object. These positions indicate where motifs are 
located within promoter sequences of a background dataset.
}
\details{
Methods should be implemented for specific object classes that store motif 
scanning results.


The `ps_hits_pos_bg` function is particularly relevant for background 
datasets, which include motif scan results computed across all promoter 
sequences.

In a PSMatrixList object, the `ps_hits_pos_bg` slot is populated only for 
background matrices (i.e., matrices derived from all promoters). 
This ensures that the pscan() function can access precomputed metrics 
without recomputing them, optimizing efficiency.
}
\examples{
full_pfm1_path <- system.file("extdata", 
                              "full_pfm1.rds", 
                              package = "PscanR")
full_pfm1 <- readRDS(full_pfm1_path)
ps_hits_pos_bg(full_pfm1)

}
