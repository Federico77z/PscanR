% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-PSMatrix.R
\name{ps_hits_pos,PSMatrix-method}
\alias{ps_hits_pos,PSMatrix-method}
\title{Get Motif Hit Positions}
\usage{
\S4method{ps_hits_pos}{PSMatrix}(x, pos_shift = 0L, withDimnames = TRUE)
}
\arguments{
\item{x}{A `PSMatrix` object.}

\item{pos_shift}{Integer. Specifies the amount to shift the position. 
Default is set to `0`.}

\item{withDimnames}{Logical, whether to include dimension names in the 
output, if they exist in the object.
Default set to `TRUE`.}
}
\value{
A named integer vector where names correspond to promoter sequence 
   identifiers, and values represent the adjusted motif hit positions.
}
\description{
Retrieves the positions of hits stored in a `PSMatrix` object. These 
positions indicate where the motifs are located in each promoter sequence. 
The positions can be shifted by a specified value to find the correspective
position in respect to the TSS.
}
\examples{
pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
pfm1 <- readRDS(pfm1_path)
ps_hits_pos(pfm1)

ps_hits_pos(pfm1, pos_shift = -200) # in respect to the TSS

}
