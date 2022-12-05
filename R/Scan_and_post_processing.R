#' Executes the Pscan algorithm on a set of sequences.
#' @param x 
#' The set of regulatory sequences from co-regulated or co-expressed genes (i.e.
#' a set of gene promoters). This must be a DNAStringSet object (see Biostrings package)
#' @param pfms 
#' An object of PSMatrixList class containing PWMs and background values. See 
#' ps_build_bg, ps_build_bg_from_file, ps_build_bg_from_table for how to create
#' PSMatrixList objects.
#' @param BPPARAM 
#' The BPPARAM used by bplapply. See BiocParallel.
#' @param BPOPTIONS 
#' The BPOPTIONS used by bplapply. See BiocParallel.
#' @return
#' Returns pfms with foreground values computed on x.
#' @export
#' @examples
#' 
pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms,type = 4)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- BiocParallel::bplapply(pfms, FUN = ps_scan, x, BG = FALSE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  
  BiocGenerics::do.call(PSMatrixList, pfms)
}

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

ps_z_table <- function(pfms)
{
  .ps_checks2(pfms)
  
  tbl <- lapply(pfms, ps_hits_z)
  
  as.matrix(as.data.frame(tbl, col.names = name(pfms)))

}