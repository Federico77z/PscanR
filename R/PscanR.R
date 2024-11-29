#' PscanR: An R Implementation for the `PScan` Algorithm
#' 
#' @description
#' `PscanR` provides an R implementation of the `Pscan` algorithm for 
#' transcription factor binding site motif analysis. 
#' It includes function to build the background for statistical analyses, 
#' to scan the promoter sequences, and data visualization.  
#' The package is designed to work with Bioconductor data objects.
#' 
#' @section Main functions:
#' \itemize{
#'    \item \code{\link{ps_build_bg_from_file}}: Build background matrices 
#'    from a file. 
#'    \item \code{\link{pscan}}: Pscan algorithm utilization on a set 
#'    of gen promoters.
#'    \item \code{\link{ps_result_table}}: Data visualization. 
#' }
#' 
#' @docType package
#' @name PscanR
#' 
#' @references 
#' `Pscan` platform: \url{http://159.149.160.88/pscan/}
#' 
#' Zambelli F, Pesole G, Pavesi G. Pscan: finding over-represented transcription 
#' factor binding site motifs in sequences from co-regulated or co-expressed 
#' genes. Nucleic Acids Res. 2009 Jul;37(Web Server issue):W247-52. 
#' doi: 10.1093/nar/gkp464. Epub 2009 May 31. PMID: 19487240; PMCID: PMC2703934.
#' 
#' @author 
#' Federico Zambelli [aut], Giulio Pavesi [cre]
#' 
#' Maintainer: Federico Zambelli <federico.zambelli@unimi.it>
NULL