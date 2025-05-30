% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-PSMatrix.R
\name{ps_hits_oligo_bg,PSMatrix-method}
\alias{ps_hits_oligo_bg,PSMatrix-method}
\title{Get Matched Oligonucleotide Sequences of a Background Dataset}
\usage{
\S4method{ps_hits_oligo_bg}{PSMatrix}(x, withDimnames = TRUE)
}
\arguments{
\item{x}{A `PSMatrix` object, typically the result of the `Pscan` algorithm.}

\item{withDimnames}{Logical, whether to include dimension names in the 
output, if they exist in the object.
Default set to `TRUE`.}
}
\value{
A character vector containing the sequences of motif matches 
    (oligonucleotides) in the banckground set.
}
\description{
Retrieves the oligonucleotide sequences that match the motif in the 
background sequences, returning a character vector.
}
\details{
This method is specifically designed for background datasets,
where motif hit oligonucleotide sequences are precomputed for all 
promoter sequences. 
The function extracts the character values from the ps_hits_oligo_bg slot 
and, if applicable, assigns sequence names using .ps_bg_seq_names().

These background metrics are particularly useful in motif enrichment 
analyses, as they allow \code{pscan()} to compare promoter 
sequences against a reference distribution without the need for 
recomputation.
}
\examples{
full_pfm1_path <- system.file("extdata", 
                              "full_pfm1.rds", 
                              package = "PscanR")
full_pfm1 <- readRDS(full_pfm1_path)
ps_hits_oligo_bg(full_pfm1)

}
