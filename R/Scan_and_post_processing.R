#' Executes the Pscan algorithm on a set of regulatory sequences.
#' 
#' This function computes the alignment scores between regulatory sequences 
#' (a set of gene promoters (`x` parameter)) and position weight matrices, which 
#' quantify potential binding affinities of transcription factors in that region. 
#' This is done throughout the Pscan algorithm. 
#' 
#' @param x A `DNASetString` object containing the set of regulatory sequences 
#'    from co-regulated or co-expressed genes (i.e. a set of gene promoters). 
#'    See the Biostrings package for details.
#'    
#' @param pfms An object of `PSMatrixList` class containing PWMs and background 
#'    values. See \code{\link{ps_build_bg}}, \code{\link{ps_retrieve_bg_from_file}}, 
#'    \code{\link{ps_build_bg_from_table}} for how to create `PSMatrixList` 
#'    objects that contain background statistics.
#' 
#' @param BPPARAM The BPPARAM used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`.
#'    If BPPARAM is not explicitly set, the default value (bpparam()) will be 
#'    used, which automatically chooses a sensible parallelization method based 
#'    on the user's system. 
#'    You can specify BPPARAM = BiocParallel::SnowParam(8) on all operating 
#'    systems, or BPPARAM = BiocParallel::MulticoreParam(8) on Unix-like
#'    systems to use, for example, 8 cores. 
#'    
#' @param BPOPTIONS The BPOPTIONS used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`. 
#'    The default is `bpoptions()`.
#'    Some useful tasks: bpoptions(progressbar = TRUE, log = TRUE). 
#'    progressbar = TRUE enables a progress bar that can be useful when 
#'    processing many tasks. log = TRUE enable logging to debug each step of
#'    the parallel tasks. 
#'    
#' @details
#' The `pscan` function performs sequence scanning using the `ps_scan` method for 
#' individual PWMs, accounting for both the forward and reverse complement strands,
#' ensuring no potential binding sites are missed.
#'    
#' @return
#' A `PSMatrixList` object in which the foreground value (the alignment scores)
#' have been computed for each sequence in `x` based on the position weight 
#' matrices in `pfms`.
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
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the PScan algorithm
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package). 
#' 
#' ps_results_table(results)
#' 
pscan <- function(x, pfms, BPPARAM=bpparam(), BPOPTIONS = bpoptions())
{
  .ps_checks(x, pfms,type = 4)
  
  x <- BiocGenerics::unique(x)
  x <- .clean_sequence(x)
  
  pfms <- BiocParallel::bplapply(
    pfms, 
    FUN = ps_scan, 
    x, 
    BG = FALSE, 
    BPPARAM=BPPARAM, 
    BPOPTIONS = BPOPTIONS)
  
  BiocGenerics::do.call(PSMatrixList, pfms)
}

#' Create a Summary Table of PscanR Results  
#' 
#' This function generates a table summarizing the Pscan results stored in a 
#' `PSMatrixList` object. It retrieves pre-computed background and foreground 
#' statistics (the alignment scores between regulatory sequences 
#' and position weight matrices) for each PWM, calculates adjusted p-values 
#' (FDR), and returns a table ordered by decreasing `ZSCORE` and 
#' increasing `P.VALUE`.  
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` function. 
#'
#' @return
#' A data.frame with matrices ordered by increasing P.VALUE and decreasing 
#' ZSCORE, containing the following columns: 
#' \itemize{
#'   \item "NAME": The names of the matrices.
#'   \item "BG_AVG": The average background score for each PWM.
#'   \item "BG_STDEV": The standard deviation of the background scores for each 
#'   PWM.
#'   \item "FG_AVG": The average foreground score for each PWM.
#'   \item "ZSCORE": The Z-score for each PWM.
#'   \item "P.VALUE": The p-value for each PWM.
#'   \item "FDR": The adjusted p-value (False Discovery Rate) calculated using 
#'   the Benjamini-Hochberg method.
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
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(12))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' table <- ps_results_table(results)
#' 
#' @seealso \code{\link{ps_generics}}
#' 
#' @export
#' @importFrom stats p.adjust
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

#' Generate a Z-score Table for Motifs
#' 
#' This function calculates the Z-score for each motif in a `PSMatrixList` 
#' object and organizes the results into a matrix. The Z-Score represents
#' the statistical significance of the alignment scores for regulatory 
#' sequences relative to background expectation.
#'
#' @param pfms A `PSMatrixList` object containing multiple Position Weight 
#'    Matrices and associated metadata (foreground and background statistics). 
#'    This object is the output of `pscan()` function. 
#'
#' @return
#' A matrix in which each column corresponds to a motif in the `PSMatrixList`
#' object, and each row to the sequence identifiers. The cell values in the 
#' matrix are Z-scores, indicating the statistical significance of the alignment 
#' between each sequence and the motif.
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
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' ps_z_table(results)
#' 
#' @export
ps_z_table <- function(pfms)
{
  .ps_checks2(pfms)
  
  tbl <- lapply(pfms, ps_hits_z)
  
 # as.matrix(as.data.frame(tbl, col.names = name(pfms)))
  
  as.matrix(as.data.frame(tbl, col.names = ID(pfms)))

}

#' Pscan Score Correlation Heatmap
#' 
#' This function generates a heatmap that visualizes the correlation of 
#' Z-scores for transcription factors (TFs) based on a specified false discovery 
#' rate threshold. It allows customization of the heatmap's appearance.  
#' 
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` function. 
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to include in the analysis. The default is `0.01`.
#' @param ... Additional user defined arguments to customize the heatmap settings, 
#'    such as color palettes or clustering object.
#'   
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts the result table and the z-score table from the pscan 
#'      algorithm result.
#'   \item Filters the result table based on the given FDR threshold.
#'   \item Uses the filtered result table to create the z-table for the 
#'      selected TFs.
#'   \item Generates the heatmap using the `pheatmap` function, 
#'      with costumizable settings.}
#' 
#'  Default settings, that can be changed by the users, are: 
#'  \itemize{
#'    \item `cluster_rows` and `cluster_cols` set to `TRUE`.
#'    \item `color` uses a blue-white-red palette.
#'    \item The main `Pscan Score Correlation Heatmap`.
#'    \item TF names are showed as column labels.
#'    }
#'
#' @return A heatmap plot showing Z-score correlations for selected 
#'    transcription factors. 
#' 
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:50]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' ps_score_correlation_map(results, FDR = 0.05)
#' 
#' @export
#' @import pheatmap 
#' @importFrom utils modifyList
ps_score_correlation_map <- function(pfms, FDR = 0.01, ...)
{
  res_table <- ps_results_table(pfms)
  z_table <- ps_z_table(pfms)
  topn <- which(res_table$FDR <= FDR)
  
  defaults <- list(cluster_rows = TRUE, 
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Pscan Score Correlation Heatmap", scale = "row", 
  show_rownames = FALSE, 
  labels_col = res_table$NAME[topn], 
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete")
  
  user_args <- list(...)
  
  final_args <- modifyList(defaults, user_args)

  tf_to_plot <- rownames(res_table)[topn]
  
  z_table_reduced <- z_table[,tf_to_plot]
  
  do.call(pheatmap, c(list(z_table_reduced), final_args))
  
  invisible(z_table_reduced)
}

#' Pscan Hits Position Heatmap
#' 
#' This function creates a heatmap visualizing the positional distribution
#' of hits based on a specified false discovery rate threshold. 
#' 
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` function. 
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to be included in the analysis. The default is set to `0.01`.
#' @param shift Numeric. A value to shift the positions of hits. 
#'    Default is set to `0`.
#' @param ... Additional user defined arguments that can be passed to 
#'    the function (e.g., the color palette) to change the default settings.
#'   
#' @details
#' The function performs the following steps: 
#' \itemize{
#'   \item Extracts the result table from the pscan algorithm result.
#'   \item Filters the result table based on the given FDR threshold.
#'   \item Creates a positional hits matrix.
#'   \item Generates the heatmap using the `pheatmap` function, 
#'      with costumizable settings.}
#' 
#' Default settings, that can be changed by the users, are:
#' \itemize{
#'   \item `cluster_rows` and `cluster_cols` set to `TRUE`.
#'   \item `color` uses a white-yellow-red palette.
#'   \item The main `Pscan Hits Position Heatmap`.
#'   }
#' 
#' @return A heatmap plot showing positional hits distribution for selected 
#'    transcription factors.
#' 
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:50]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' ps_hitpos_map(results)
#' 
#' 
#' @export
#' @import pheatmap 
#' @importFrom grDevices colorRampPalette
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
  
  pos_mat <- matrix(data = NA, nrow = ps_fg_size(results[[1]]), 
                    ncol = length(topn))
  
  for(v in topn)
    pos_mat[,v] <- ps_hits_pos(results[[row.names(res_table)[v]]], 
                               pos_shift = shift)
  
  colnames(pos_mat) <- res_table$NAME[topn]
  rownames(pos_mat) <- ps_seq_names(pfms[[1]])
  
  do.call(pheatmap, c(list(pos_mat), final_args))
  
  invisible(pos_mat)
}

#' Pscan Density Plot of Hits along Promoters 
#' 
#' This function creates a density plot representing the distribution of hits 
#' along the promoter region based on their position and score. 
#' 
#' @param pfm A Position Frequency Matrix, result of the Pscan algorithm. 
#'
#' @param shift Numeric value specifying the positional shift applied to the hit 
#'    positions. Default is `0`.
#' @param st Score threshold used to filter hits. Can be a numeric value to set 
#'    the threshold directly, or a character:
#'    \itemize{
#'      \item `all`: the threshold is set to `0` (so, no threshold. All the hits
#'      are evaluated).
#'      \item `loose`: uses the background average as threshold.
#'      \item `strict`: uses the background average plus the background standard 
#'      deviation as threshold.}
#'    Default is set to loose.
#'      
#' @return A density plot showing the distribution of hits along the promoter 
#'    region.
#'      
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:50]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' pfm1 <- results[[1]]
#' ps_density_plot(pfm1)
#' 
#' @export
#' @importFrom grDevices rgb
#' @importFrom graphics abline
#' @importFrom graphics polygon
#' @importFrom graphics text
#' @importFrom stats density
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
  text(peak, max(density_hits$y), 
       labels = paste("\tMode:", round(peak)), pos = 4)
  
}

#' Bubble chart of Pscan Hits Score vs Position
#' 
#' This function visualizes the Bubble Chart to visualize the relationship
#' between position and score of the identified sites in a `PSMatrix` object. 
#' 
#' @param pfm An object of class `PSMatrix`. It must be processed by the PscanR 
#'    algorithm.
#'    
#' @return A ggplot object representing a bubble chart. The x-axis corresponds 
#' to the hits position, the y-axis corresponds to the hits score, and the 
#' size of the bubbles represents the count of occurrences. 
#' 
#' @details
#' The function computes the count of occurrences for each combination of 
#' PS hits score and position from the input `PSMatrix` and visualizes 
#' this data using a bubble chart. The user can modify the color of bubbles. 
#' Default is `blue`.
#' 
#'
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:50]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' pfm1 <- results[[1]]
#' ps_score_position_BubbleChart(pfm1)
#' 
#' @export
#' @import dplyr 
#' @import ggplot2
ps_score_position_BubbleChart <- function(pfm, bubble_color = 'blue')
{
  data <- ps_hits_table(pfm)
  
  data_sum <- data %>% 
    group_by(!!sym(colnames(data)[1]), !!sym(colnames(data)[2])) %>%
    summarise(Count = n(), .groups = "drop")
  
  ggplot(data_sum, aes(x = !!sym(colnames(data_sum)[2]), 
                       y = !!sym(colnames(data_sum)[1]), size = Count)) +
    geom_point(alpha=0.5, color = bubble_color) +
    scale_size_continuous(breaks = sort(unique(data_sum$Count)), 
                          guide = guide_legend(title = "Occurrences")) +
    labs(x = "PS Hits Position", y = "PS Hits Score", 
         title = paste(pfm@name,"Bubble Chart of Score vs Position Hits")) +
    theme_minimal()
  
}

#' Density Plot of Distances between Identified Sites in two PSMatrix Object
#' 
#' This function visualizes the density plot of distances between identified 
#' sites in two `PSMatrix` object. It allows to filter the identified sites
#' based on a specified threshold value. 
#' 
#'
#' @param M1 An object of class `PSMatrix`. It must be processed by the PscanR 
#'    algorithm. 
#' @param M2 An object of class `PSMatrix`. It must be processed by the PscanR 
#'    algorithm. 
#' @param st1 Score threshold used to filter hits for the first PSMatrix . 
#'    Can be a numeric value to set the threshold directly, or a character:
#'    \itemize{
#'      \item `all`: the threshold is set to `0` (so, no threshold. All the hits
#'      are evaluated).
#'      \item `loose`: uses the background average as threshold.
#'      \item `strict`: uses the background average plus the background standard 
#'      deviation as threshold.}
#'    Default is set to loose.
#' @param st2 Score threshold used to filter hits for the second PSMatrix. 
#'    Can be a numeric value to set the threshold directly, or a character, as 
#'    for st1. 
#'    Default is set to loose. 
#'
#' @return A density plot visualizing the distances between identified sites 
#'    in \code{M1} and \code{M2}.
#' 
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_bg_std_dev}}, 
#' \code{\link{ps_hits_score}}, \code{\link{ps_hits_pos}}
#' 
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:50]
#' 
#' # Load JASPAR motif matrices for vertebrates
#' J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
#' load(J2020_path)
#' 
#' bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
#'                        package = "PscanR")
#' J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' pfm1 <- results[[1]]
#' pfm2 <- results[[2]]
#' ps_density_distances_plot(pfm1, pfm2, 'all', 'loose')
#' 
ps_density_distances_plot <- function(M1, M2, st1 = ps_bg_avg(M1), st2 = ps_bg_avg(M2))
{
  if(class(M1) != "PSMatrix" || class(M2) != "PSMatrix") {
    stop("Both object must be of class PSMatrix")
  }
  
  if (is.character(st1)){
    st1 <- switch(st1,
                  "all" = 0,
                  "loose" = ps_bg_avg(M1),
                  "strict" = ps_bg_avg(M1) + ps_bg_std_dev(M1),
                  {
                    warning("Invalid value for st1, reverting to loose")
                    ps_bg_avg(M1)
                  })
  }
  
  if (is.character(st2)) {
    st2 <- switch(st2,
                  "all" = 0,
                  "loose" = ps_bg_avg(M2),
                  "strict" = ps_bg_avg(M2) + ps_bg_std_dev(M2),
                  {
                    warning("Invalid value for st2, reverting to loose")
                    ps_bg_avg(M2)
                  })
  }
  
  scores1 <- ps_hits_score(M1)
  g_scores1 <- scores1 >= st1
  sum_g1 <- sum(g_scores1)
  
  scores2 <- ps_hits_score(M2)
  g_scores2 <- scores2 >= st2
  sum_g2 <- sum(g_scores2)
  
  hits1 <- ps_hits_pos(M1, pos_shift = ncol(M1)/2)[g_scores1]
  hits2 <- ps_hits_pos(M2, pos_shift = ncol(M2)/2)[g_scores2]
  
  hits1 <- hits1[names(hits1) %in% names(hits2)]
  hits2 <- hits2[names(hits2) %in% names(hits1)]
  
  distances <- hits1 - hits2
  density_distances <- density(distances)
  
  plot(density_distances, 
       main = paste(M1@name, 'and', M2@name, 'distances density plot'),
       xlab = "Distances between the identified sites",
       ylab = "Density",
       col = "blue",
       lwd = 2)
  polygon(density_distances, col = rgb(0,0,1,0.1), border = NA)
  
  #peak <- density_distances$x[which.max(density_distances$y)]
  #abline(v = peak, col = "gray", lty = 2, lwd = 2)
  #text(peak, max(density_distances$y), 
  #     labels = paste("\tMode:", round(peak)), pos = 4)
}


#' Generate PSMatrixList from Background and JASPAR Matrix
#' 
#' This function generates an object of class `PSMatrixList` using a specified 
#' JASPAR matrix, an organism identifier, and a promoter region range. 
#' It is designed to work with only JASPAR database and UCSC annotations. 
#'
#' @param JASPAR_matrix A character string specifying the JASPAR database 
#'    version. You can choose between 'JASPAR2020', 'JASPAR2022', or 
#'    'JASPAR2024'
#' @param org The accepted values are: `hs` (Homo sapiens), `mm` (Mus musculus), 
#'    `at` (Arabidopsis thaliana), `sc` (Saccharomyces cerevisiae), `dm`
#'    (Drosophila melanogaster). 
#' @param prom_reg A numeric vector of two integers indicating the range 
#'    of the promoter region. You can choose between:
#'    \itemize{
#'      \item 200 base pairs upstream and 50 downstream base pair from the TSS 
#'      (c(-200, 50))
#'      \item 450 bp upstream and 50 downstream (c(-450, 50))
#'      \item 500 bp upstream and 0 downstream (c(-500, 0))
#'      \item 950 bp upstream and 50 downstream (c(-950, 50))
#'      \item 1000 bp upstream and 0 downstream (c(-1000, 0))
#'      } 
#'
#' @return A `PSMatrixList` object created from the specified background 
#' file and JASPAR matrix
#' 
#' @export
#'
#' @examples
#' generate_psmatrixlist_from_background('Jaspar2020', 'hs', c(-200,50))
#' 
generate_psmatrixlist_from_background <- function(JASPAR_matrix, org, prom_reg, assembly){
  
  # Rivedi funzione per l'assembly
  
  organism_map <- list(
    "hs" = assembly,
    "mm" = assembly,
    "at" = "at10",
    "sc" = "sacCer3",
    "dm" = "dm6"
  )
  
  if(!(org %in% names(organism_map))){
    stop("Invalid organism acronym. Choose between: 'hs', 'mm', 'at', 'sc', 'dm'.")
  }
  
  org <- organism_map[[org]]
  
  J_name <- toupper(JASPAR_matrix)
  
  if(!(J_name == 'JASPAR2020' || J_name == 'JASPAR2022' || J_name == 'JASPAR2024')){
    stop("Invalid database. Choose between: 'JASPAR2020', 'JASPAR2022', 'JASPAR2024'
         (non case sensitive).")
  }
  
  version <- substr(JASPAR_matrix,7, 10)
  
  if(length(prom_reg) != 2 || !is.numeric(prom_reg)) {
    stop("Promoter region must be a numeric vector with two integer.")
  }
  p_up <- abs(prom_reg[1])
  p_down <- abs(prom_reg[2])
  
  # Rivedi funzione per gli altri org, il nome dei file sarà diverso!
  
  file_name <- paste0('J', version, '_', org, '_', p_up, 'u_', p_down, 'd_UCSC.psbg.txt')

  
  BG_path <- system.file("extdata/BG_scripts", file_name, package = 'PscanR')
  
  opts <- list()
  opts[["collection"]] <- "CORE" 
  opts[["tax_group"]] <- "vertebrates" 
  
  # ATTENZIONE: non ancora testata per JASPAR2024
  
  if(J_name == 'JASPAR2020'){
    J_matrix <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) 
  }
  if(J_name == 'JASPAR2022'){
    J_matrix <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
  }
  if(J_name == 'JASPAR2024'){
    httr::set_config(config(ssl_verifypeer = 0L))
    JASPAR2024 <- JASPAR2024()
    JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
    J_matrix <- TFBSTools::getMatrixSet(JASPARConnect, opts)
  }
  
  J_PSBG <- ps_retrieve_bg_from_file(BG_path, J_matrix)
  return(J_PSBG)
}