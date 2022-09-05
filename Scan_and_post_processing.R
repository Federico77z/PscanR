pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms,type = 4)
  
  x <- BiocGenerics::unique(x)
  
  pfms <- bplapply(pfms, FUN = ps_scan, x, BG = FALSE, BPPARAM=BPPARAM, BPOPTIONS = BPOPTIONS)
  
  do.call(PSMatrixList, pfms)
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