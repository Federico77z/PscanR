#' PscanR: An R Implementation for the `Pscan` Algorithm
#'
#' @description
#' `PscanR` provides an R implementation of the `Pscan` algorithm for
#' transcription factor binding site motif analysis.
#' The package is designed to work with Bioconductor data objects:
#' it extends PFMatrix and PFMatrixList from
#' the TFBSTools package for transcription factor binding profiles and
#' employs Biostrings DNAStringSet for promoter sequences.
#' PscanR supports parallel execution through the BiocParallel package.
#' PscanR includes functions to build the background for a set of promoter
#' sequences, to scan a foreground set of promoter sequences, and to visualize
#' and plot the results.
#'
#' @section Main functions:
#' \itemize{
#'    \item \code{\link{generate_psmatrixlist_from_background}}:
#'    Build background matrices from pre-computed files.
#'    \item \code{\link{pscan}}:
#'    Pscan algorithm utilization on a set of gene promoters.
#'    \item \code{\link{ps_results_table}}:
#'    Data visualization.
#' }
#'
#' @examples
#' motifs <- readRDS(system.file(
#'     "extdata", "J2020.rds", package = "PscanR"
#' ))
#' background <- ps_retrieve_bg_from_file(system.file(
#'     "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'     package = "PscanR"
#' ), motifs)
#' # NRF1 and CTCF illustrate enriched and background-like motif scores.
#' background <- background[c("MA0506.1", "MA0139.1")]
#' # Historical NRF1-example promoters; the original experiment is unknown.
#' sequences <- readRDS(system.file(
#'     "extdata", "prom_seq.rds", package = "PscanR"
#' ))
#' # Adjusted p-values cover only these two motifs.
#' ps_results_table(pscan(sequences, background))
#'
#' @docType package
#' @name PscanR
#'
#' @references
#' `Pscan` Web: \url{http://www.beaconlab.it/pscan/}
#'
#' Zambelli F, Pesole G, Pavesi G. Pscan: finding over-represented
#' transcription factor binding site motifs in sequences from co-regulated or
#' co-expressed genes.
#' Nucleic Acids Res. 2009 Jul;37(Web Server issue):W247-52.
#' doi: 10.1093/nar/gkp464. Epub 2009 May 31. PMID: 19487240;
#' PMCID: PMC2703934.
#'
#' @author
#' Federico Zambelli [aut, cre], Giulio Pavesi [aut]
#'
#' Maintainer: Federico Zambelli <federico.zambelli@unimi.it>
#'
#' @importFrom methods as callNextMethod is new validObject
#' @importFrom stats pnorm sd setNames
#' @importFrom utils read.table write.table
#' @importFrom BiocParallel bplapply bpoptions
#' @importFrom Biostrings DNAStringSet
#' @keywords internal
"_PACKAGE"

NULL
