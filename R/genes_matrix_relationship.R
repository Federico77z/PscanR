#' Filter and Scan Promoter Regions with Pscan
#' 
#' This function filters input promoter sequences based on their affinity 
#' with a specified JASPAR matrix before applying the Pscan algorithm with the 
#' background reference.
#' 
#' @param prom_seq A DNAStringSet object. The promoter sequences to be filtered.
#' @param Jmatrix A PSMatrix object. A JASPAR matrix containing background 
#'    average and standard deviation statistics, used to filter promoter 
#'    sequences.
#' @param n A numeric value determining the threshold for filtering. The 
#'    threshold is computed as n*ps_bg_std_dev(Jmatrix) + ps_bg_avg(Jmatrix). 
#'    Positive values retain sequences with higher affinity, while 
#'    negative values retain sequences with lower affinity. Default is 1.
#' @param background A PSMatrixList object. It contains PWMs and background statistics.
#'    For background statistics we refer to the standard deviation and average 
#'    of hits scores when the background (set of promoters of all the transcript 
#'    in the organism of study) is scanned with the position weight matrices. 
#'    This is used to assert the statistical enrichment of motif occurrences 
#'    in co-expressed or co-regulated genes.
#'    See \code{\link{ps_build_bg}}, \code{\link{ps_retrieve_bg_from_file}}, 
#'    \code{\link{ps_build_bg_from_table}}, \code{\link{generate_psmatrixlist_from_background}}
#'    for how to create `PSMatrixList` objects that contain background statistics. 
#' @param BPPARAM The BPPARAM used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`.
#'    If BPPARAM is not explicitly set, the default value (bpparam()) will be 
#'    used, which automatically chooses a sensible parallelization method based 
#'    on the user's system. 
#'    You can specify BPPARAM = BiocParallel::SnowParam(8) on all operating 
#'    systems, or BPPARAM = BiocParallel::MulticoreParam(8) on Unix-like
#'    systems to use, for example, 8 cores. 
#' @param BPPOPTIONS The BPOPTIONS used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`. 
#'    The default is `bpoptions()`.
#'    Some useful tasks: bpoptions(progressbar = TRUE, log = TRUE). 
#'    progressbar = TRUE enables a progress bar that can be useful when 
#'    processing many tasks. log = TRUE enable logging to debug each step of
#'    the parallel tasks.  description
#'    
#' @details
#'    This function:
#'    \itemize{
#'       \item scans promoter sequences using the `Jmatrix` PWM.
#'       \item retains only those that meet the threshold criterion. 
#'       \item scans the filtered sequences across the background PWMs.
#'    }
#' 
#' @seealso \code{\link{pscan}}
#' 
#' @return A `PSMatrixList` object in which the foreground values 
#' (the alignment scores) have been computed for each sequence in `prom_seq`, 
#' that has passed the threshold, based on the position weight matrices 
#' in `background`. If no sequences pass the filter, a warning is issued, and  
#' the function returns `NULL`.  
#' @export
PscanFiltered <- function(prom_seq, Jmatrix, n = 1, background, 
                          BPPARAM=bpparam(), BPPOPTIONS = bpoptions()){
  
  if (!is(prom_seq, "DNAStringSet"))
    stop("Invalid input: 'prom_seq' must be a DNAStringSet")
  if (!is(Jmatrix, "PSMatrix"))
    stop("Invalid input: 'Jmatrix' must be a PSMatrix")
  if (!is(background, "PSMatrixList"))
    stop("Invalid input: 'background' must be a PSMatrixList")
  
  threshold <- ps_bg_avg(Jmatrix) + n * ps_bg_std_dev(Jmatrix)
  
  rc_matrix <- reverseComplement(Jmatrix) 
  
  Margs <- list(numx = as.numeric(Matrix(Jmatrix)), 
                numx_rc = as.numeric(Matrix(rc_matrix)),
                ncolx = (0:(ncol(Matrix(Jmatrix)) - 1))*length(.PS_ALPHABET(Jmatrix)), 
                AB = .PS_ALPHABET(Jmatrix))
  
  res <- mapply(.ps_scan_s, list(Jmatrix), as.character(prom_seq), MoreArgs = Margs)
  
  Score <- as.numeric(res["score",])
  
  filtered_prom_seq <- character()
  if(n > 0)
    filtered_prom_seq <- prom_seq[Score >= threshold]
  if(n < 0)
    filtered_prom_seq <- prom_seq[Score <= threshold]
  
  if (length(filtered_prom_seq) == 0) {
    warning("No sequence satisfy the filter criterium")
    return(NULL)
  }
  
  pfms <- BiocParallel::bplapply(
    background, 
    FUN = ps_scan, 
    filtered_prom_seq, 
    BG = FALSE, 
    BPPARAM=BPPARAM, 
    BPPOPTIONS = BPPOPTIONS)
  
  BiocGenerics::do.call(PSMatrixList, pfms)
  
}
