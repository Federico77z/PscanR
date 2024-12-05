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
#' 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR")
#' J2020_PSBG <- ps_build_bg_from_file(bg_path, J2020)
#' 
#' # Execute the PScan algorithm
#' results <- pscan(prom_seq, J2020_PSBG, BPPARAM = BiocParallel::MulticoreParam(1))
#' results
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
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR")
#' J2020_PSBG <- ps_build_bg_from_file(bg_path, J2020)
#' 
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, BPPARAM = BiocParallel::MulticoreParam(1))
#' # Use `SnowParam` on windows 
#' table <- ps_results_table(results)
#' 
#' @seealso \code{\link{ps_generics}}
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
  fdr_v <- p.adjust(pv_v, method = "BH")

  tbl <- data.frame("NAME" = name(pfms), "BG_AVG" = bg_v, "BG_STDEV" = std_v, 
             "FG_AVG" = fg_v, "ZSCORE" = zs_v, 
             "P.VALUE" = pv_v, "FDR" = fdr_v, row.names = ID(pfms))
  
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
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR")
#' J2020_PSBG <- ps_build_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, BPPARAM = BiocParallel::MulticoreParam(1))
#' # Use `SnowParam` on windows
#' z_score <- ps_z_table(results)
#' 
#' @export
ps_z_table <- function(pfms)
{
  .ps_checks2(pfms)
  
  tbl <- lapply(pfms, ps_hits_z)
  
 # as.matrix(as.data.frame(tbl, col.names = name(pfms)))
  
  as.matrix(as.data.frame(tbl, col.names = ID(pfms)))

}

#' @export
#' @import pheatmap 

ps_score_correlation_map <- function(pfms, FDR = 0.01, ...)
{
  res_table <- ps_results_table(pfms)
  z_table <- ps_z_table(pfms)
  topn <- which(res_table$FDR <= FDR)
  
  defaults <- list(cluster_rows = TRUE, 
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Pscan Score Correlation Heatmap", scale = "row", show_rownames = FALSE, 
  labels_col = res_table$NAME[topn], 
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete")
  
  user_args <- list(...)
  
  final_args <- modifyList(defaults, user_args)

  tf_to_plot <- rownames(res_table)[topn]
  
  z_table_reduced <- z_table[,tf_to_plot]
  
  do.call(pheatmap, c(list(z_table_reduced), final_args))
  
}

#' @export
#' @import pheatmap 
ps_hitpos_map <- function(pfms, FDR = 0.01, shift = 0, ...)
{
  res_table <- ps_results_table(pfms)
  
  topn <- which(res_table$FDR <= FDR)
  
  defaults <- list(cluster_rows = TRUE,  
                   cluster_cols = TRUE,  
                   color = colorRampPalette(c("white", "yellow", "red"))(100),  
                   main = "Pscan Hits Position Heatmap",
                   fontsize = 10, show_rownames = FALSE, scale = "none",
                   clustering_distance_rows = "manhattan",
                   clustering_distance_cols = "manhattan", 
                   clustering_method = "average"
  )
  
  user_args <- list(...)
  
  final_args <- modifyList(defaults, user_args)
  
  pos_mat <- matrix(data = NA, nrow = ps_fg_size(results[[1]]), ncol = length(topn))
  
  for(v in topn)
    pos_mat[,v] <- ps_hits_pos(results[[row.names(res_table)[v]]], pos_shift = shift)
  
  colnames(pos_mat) <- res_table$NAME[topn]
  rownames(pos_mat) <- ps_seq_names(pfms[[1]])
  
  do.call(pheatmap, c(list(pos_mat), final_args))
}

#' @export
ps_density_plot <- function(pfm, shift = 0, st = ps_bg_avg(pfm))
{
  # st = score threshold. It can be passed as a numeric value
  # or as one of three characters "all", "loose", "strict".
  
  if(is.character(st))
  {
    if(st == "all")
      st <- 0
    else if (st == "loose")
      st <- ps_bg_avg(pfm)
    else if (st == "strict")
      st <- ps_bg_avg(pfm) + ps_bg_std_dev(pfm)
    else {
      warning("Invalid value for st, reverting to loose")
      st <- ps_bg_avg(pfm)
    }
  }
  
  scores <- ps_hits_score(pfm)
  g_scores <- scores >= st
  sum_g <- sum(g_scores)
  
  density_hits <- density(ps_hits_pos(pfm, pos_shift = shift)[g_scores])
  
  plot(density_hits, 
       main = paste(name(pfm), "hits density on", sum_g, "promoters"),
       xlab = "Position along promoters",
       ylab = "Density",
       col = "blue",
       lwd = 2)
  
  polygon(density_hits, col = rgb(0, 0, 1, 0.1), border = NA)
  
  peak <- density_hits$x[which.max(density_hits$y)]
  abline(v = peak, col = "gray", lty = 2, lwd = 2)
  text(peak, max(density_hits$y), labels = paste("\tMode:", round(peak)), pos = 4)
  
}