% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{ps_bg_seq_names}
\alias{ps_bg_seq_names}
\title{Retrieve Sequence Names for the Background Dataset (Generic Function)}
\usage{
ps_bg_seq_names(x, ...)
}
\arguments{
\item{x}{An object containing sequence information.}

\item{...}{Additional arguments passed to specific methods.}
}
\value{
\itemize{
  \item If `x` is a `PSMatrix` object: a **character vector** of sequence 
  names corresponding to the analyzed promoter regions.
  \item If `x` belongs to another supported class, the return format may 
  vary.
}
}
\description{
`ps_bg_seq_names` is a **generic function** that extracts the sequence names 
or identifiers from an object containing all the promoter sequence data for 
a specific organism.
}
\details{
Methods should be implemented for specific object classes storing sequence 
data.


The `ps_bg_seq_names` function is particularly relevant for background 
datasets, which include motif scan results computed across all promoter 
sequences.

In a PSMatrixList object, the `ps_bg_seq_names` slot is populated only for 
background matrices (i.e., matrices derived from all promoters). 
This ensures that the pscan() function can access precomputed metrics 
without recomputing them, optimizing efficiency.
}
\examples{
full_pfm1_path <- system.file("extdata", 
                              "full_pfm1.rds", 
                              package = "PscanR")
full_pfm1 <- readRDS(full_pfm1_path)
ps_bg_seq_names(full_pfm1)

}
