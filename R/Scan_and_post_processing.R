#' Executes the Pscan algorithm on a set of regulatory sequences.
#' 
#' This function computes the foreground values of a set of gene promoters (`X`)
#' using the PScan algorithm. Foreground values are computed by applying the 
#' PWM values stored in `pfms`. 
#' 
#' @param x A `DNASetString` object containing the set of regulatory sequences 
#'    from co-regulated or co-expressed genes (i.e. a set of gene promoters). 
#'    (see Biostrings package).
#'    
#' @param pfms An object of PSMatrixList class containing PWMs and background 
#'    values. See ps_build_bg, ps_build_bg_from_file, ps_build_bg_from_table for 
#'    how to create PSMatrixList objects.
#' 
#' @param BPPARAM The BPPARAM used by bplapply. See BiocParallel.
#'    This argument is passed to `BiocParallel::bplapply`
#'    The default is `bpparam()`.
#'    
#' @param BPOPTIONS The BPOPTIONS used by bplapply. See BiocParallel.
#'    This argument is passed to `BiocParallel::bplapply`. 
#'    The default is `bpoptions()`.
#'    
#' @return
#' A `PSMatrixList` object in which the foreground value have been computed 
#' for each sequence in `X` based on the position weight matrices in `pfms`.
#' 
#' @export
#' 
#' @examples
#' # Get the regulatory sequences 
#' target <- read.csv("../PscanR/example_file_nfkb100.txt", header = FALSE) 
#' txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")
#' seqlevels(txdb) <- seqlevels(txdb)[1:24] # Use only canonical chromosomes
#' prom_rng <- promoters(txdb, upstream = 500, downstream = 0, use.names = TRUE)
#' prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)
#' target_prom_rng <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]
#' prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, target_prom_rng) 
#' 
#' opts <- list(collection = "CORE", tax_group = "vertebrates")
#' J2020 <- getMatrixSet(JASPAR2020, opts)
#' 
#' J2020_PSBG <- ps_build_bg_from_file(
#'   "../BG_scripts/J2020_hg38_500u_0d_UCSC.psbg.txt", 
#'   J2020
#' )
#' 
#' # Execute the PScan algorithm
#' results <- pscan(prom_seq, J2020_PSBG)
#' 
#' 
pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms,type = 4)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- BiocParallel::bplapply(
    pfms, 
    FUN = ps_scan, 
    x, 
    BG = FALSE, 
    BPPARAM=BPPARAM, 
    BPOPTIONS = BPOPTIONS)
  
  BiocGenerics::do.call(PSMatrixList, pfms)
}

#' View result in table format  
#' 
#' This function generates a table of results from a `PSMatrixList` object, 
#' ordered by decreasing P.VALUE and ZSCORE. For each PWMs, it computes 
#' background and #' foreground statistics and returns a table summarizing 
#' these values. 
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` function. 
#'
#' @return
#' A data.frame with matrices ordered by decreasing P.VALUE and ZSCORE, 
#' containig the following columns: 
#' \itemize{
#'   \item "NAME": The names of the matrices.
#'   \item "BG_AVG": The average background score for each PWM.
#'   \item "BG_STDEV": The standard deviation of the background scores for each 
#'   PWM.
#'   \item "FG_AVG": The average foreground score for each PWM.
#'   \item "ZSCORE": The Z-score for each PWM.
#'   \item "P.VALUE": The p-value for each PWM.
#' }
#'
#' @examples
#' 
#' # Get the regulatory sequences 
#' target <- read.csv("../PscanR/example_file_nfkb100.txt", header = FALSE) 
#' txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")
#' seqlevels(txdb) <- seqlevels(txdb)[1:24] # Use only canonical chromosomes
#' prom_rng <- promoters(txdb, upstream = 500, downstream = 0, use.names = TRUE)
#' prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)
#' target_prom_rng <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]
#' prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, target_prom_rng) 
#' 
#' opts <- list(collection = "CORE", tax_group = "vertebrates")
#' J2020 <- getMatrixSet(JASPAR2020, opts)
#' 
#' J2020_PSBG <- ps_build_bg_from_file(
#'   "../BG_scripts/J2020_hg38_500u_0d_UCSC.psbg.txt", 
#'   J2020
#' )
#' 
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, BPPARAM = MultiCoreParam(24))
#' # Use `SnowParam` on windows 
#' table <- ps_result_table(result)
#' View(table)
#' 
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_bg_std_dev}},
#'    \code{\link{ps_fg_avg}}, \code{\link{ps_zscore}}, \code{\link{ps_pvalue}}
#' 
#' @export
ps_results_table <- function(pfms)
{
  
  .ps_checks2(pfms)
  
  bg_v <- vapply(pfms, ps_bg_avg, numeric(length = 1L))
  std_v <- vapply(pfms, ps_bg_std_dev, numeric(length = 1L))
  fg_v <- vapply(pfms, ps_fg_avg, numeric(length = 1L))
  zs_v <- vapply(pfms, ps_zscore, numeric(length = 1L))
  pv_v <- vapply(pfms, ps_pvalue, numeric(length = 1L))

  tbl <- data.frame("NAME" = name(pfms), "BG_AVG" = bg_v, "BG_STDEV" = std_v, 
             "FG_AVG" = fg_v, "ZSCORE" = zs_v, 
             "P.VALUE" = pv_v, row.names = ID(pfms))
  
  tbl[with(tbl, order(P.VALUE, ZSCORE, decreasing = c(FALSE,TRUE))),]
}

#' Generates Z-score table
#' 
#' This function calculates the z-score for each motif in a `PSMatrixList` 
#' object and organizes the result into a matrix. 
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` function. 
#'
#' @return
#' A matrix in which each column corresponds to a motif in the `PSMatrixList`
#' object, and each row to the sequence identifiers. 

#'
#' @examples
#' # Get the regulatory sequences 
#' target <- read.csv("../PscanR/example_file_nfkb100.txt", header = FALSE) 
#' txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")
#' seqlevels(txdb) <- seqlevels(txdb)[1:24] # Use only canonical chromosomes
#' prom_rng <- promoters(txdb, upstream = 500, downstream = 0, use.names = TRUE)
#' prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)
#' target_prom_rng <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]
#' prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, target_prom_rng) 
#' 
#' opts <- list(collection = "CORE", tax_group = "vertebrates")
#' J2020 <- getMatrixSet(JASPAR2020, opts)
#' 
#' J2020_PSBG <- ps_build_bg_from_file(
#'   "../BG_scripts/J2020_hg38_500u_0d_UCSC.psbg.txt", 
#'   J2020
#' )
#' 
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, BPPARAM = MultiCoreParam(24))
#' # Use `SnowParam` on windows
#' z_score <- ps_z_table(result)
#' View(z_score)
#' 
#' @export
ps_z_table <- function(pfms)
{
  .ps_checks2(pfms)
  
  tbl <- lapply(pfms, ps_hits_z)
  
 # as.matrix(as.data.frame(tbl, col.names = name(pfms)))
  
  as.matrix(as.data.frame(tbl, col.names = ID(pfms)))

}
