#' PscanR: An R Implementation for the `Pscan` Algorithm
#' 
#' @description
#' `PscanR` provides an R implementation of the `Pscan` algorithm for 
#' transcription factor binding site motif analysis. 
#' The package is designed to work with Bioconductor data objects: 
#' it extends PFMatrix and PFMatrixList from 
#' the TFBSTools package for transcription factors binding profiles and
#' employs Biostring DNAStringSet for promoter sequences. 
#' PscanR supports multithreading through the BiocParallel package.
#' PscanR includes functions to build the background for a set of promoter sequences,
#' to scan a foreground set of promoter sequences, and to visualize and plot the results.  
#' 
#' @section Main functions:
#' \itemize{
#'    \item \code{\link{ps_retrieve_bg_from_file}}: Build background matrices 
#'    from a file. 
#'    \item \code{\link{pscan}}: Pscan algorithm utilization on a set 
#'    of gen promoters.
#'    \item \code{\link{ps_results_table}}: Data visualization. 
#' }
#' 
#' @docType package
#' @name PscanR
#' 
#' @references 
#' `Pscan` Web: \url{http://159.149.160.88/pscan/}
#' 
#' Zambelli F, Pesole G, Pavesi G. Pscan: finding over-represented transcription 
#' factor binding site motifs in sequences from co-regulated or co-expressed 
#' genes. Nucleic Acids Res. 2009 Jul;37(Web Server issue):W247-52. 
#' doi: 10.1093/nar/gkp464. Epub 2009 May 31. PMID: 19487240; PMCID: PMC2703934.
#' 
#' @author 
#' Federico Zambelli [cre], Diana Betelli [aut], Giulio Pavesi [cre]
#' 
#' Maintainer: Federico Zambelli <federico.zambelli@unimi.it>
#' 
#' @importFrom methods as callNextMethod is new validObject
#' @importFrom stats sd setNames
#' @importFrom utils read.table write.table
#' @importFrom BiocParallel bplapply bpoptions bpparam
#' @importFrom BSDA z.test
#' @importFrom Biostrings DNAStringSet
NULL
