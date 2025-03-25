#' Executes the Pscan algorithm on a set of regulatory sequences.
#' 
#' This function computes the alignment scores between regulatory sequences 
#' (a set of gene promoters (`x` parameter)) and position weight matrices, which 
#' quantify potential binding affinities of transcription factors in that region.
#' These matrices can be sourced from public databases such as JASPAR.
#' The scanning is performed throughout the Pscan algorithm. 
#' 
#' @param x A `DNAStringSet` object containing the set of regulatory sequences 
#'    from co-regulated or co-expressed genes (i.e. a set of gene promoters). 
#'    See the Biostrings package for details.
#'    
#' @param pfms A `PSMatrixList` object containing PWMs and background statistics.
#'    For background statistics we refer to the standard deviation and average 
#'    of hits scores when the background (set of promoters of all the transcript 
#'    in the organism of study) is scanned with the position weight matrices. 
#'    This is used to assert the statistical enrichment of motif occurrences 
#'    in co-expressed or co-regulated genes.
#'    See \code{\link{ps_build_bg}}, \code{\link{ps_retrieve_bg_from_file}}, 
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
#' The `pscan` function performs sequence scanning using the `ps_scan` method 
#' for individual PWMs, accounting for both the forward and reverse complement 
#' strands, ensuring no potential binding sites are missed.
#' The method extract all the k-mers substring of the input sequences with the 
#' same length of the motif and evaluates a score for both the forward and 
#' reverse strand based on how well the k-mer matches the Transcription Factor
#' Binding Motif provided by the PWM. For each input sequence scanned with 
#' individual PWM, only the k-mer with the highest score is selected. 
#' The method will populate the following PSMatrix slots for each input sequence: 
#' \itemize{
#'    \item ps_hits_score: the highest matching value found for each sequence. 
#'    \item ps_hits_pos: position of the best TFBS match
#'    \item ps_hits_strand: the DNA strand of the TFBS ('+' for forward and '-'
#'    for reverse)
#'    \item ps_hits_oligo: the sequence of the binding site} 
#'    
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples. 
#'    
#' @return
#' A `PSMatrixList` object in which the foreground values (the alignment scores)
#' have been computed for each sequence in `x` based on the position weight 
#' matrices in `pfms`. 
#' 
#' @export
#' 
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS.
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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

#' Extract Pre-Computed Metrics from a PSMatrixList Object
#' 
#' This function is an optimized version of `pscan()`, designed to significantly
#' reduce computational time by extracting relevant metrics from a pre-computed
#' background PSMatrixList. 
#' 
#' If a fully computed PSMatrixList is available, containing scores, positions, 
#' oligonucleotide sequences, and strand information for all promoter sequences 
#' of a given organism (scanned using a Position Weight Matrix), this function 
#' allows users to efficiently retrieve relevant values for specific sequences 
#' of interest, without having to recalculate them with the `pscan()` function. 
#' 
#' @param ID A character vector containing the sequence identifiers (transcript 
#'    IDs) of the sequences to be scanned with the PWMs.
#' @param full_pfms The complete background PSMatrixList generated by
#'    ps_build_bg(). To preserve key metrics, ensure the object is saved using 
#'    save() and the flag fullBG is set to `TRUE`. 
#'    
#' @return A `PSMatrixList` object in which the alignment scores and related 
#'    metrics have been retrieved from `full_pfms` for each sequence in `ID`, 
#'    based on the background data.
#' 
#' @export
pscan_full_bg <- function(ID, full_pfms)
{  
   .ps_checks(x, pfms,type = 4)
   if(!is.character(ID))
     stop('ID must be a character vector containing transcript identifiers')
  
   all_seq_ID <- full_pfms@transcriptIDLegend
   
   # Remove transcript extension for any name 
   names(all_seq_ID) <- sub("\\..*$", "", names(all_seq_ID))
   all_seq_ID <- sub("\\..*$", "", all_seq_ID)
   ID <- sub("\\..*$", "", ID)
   
   # Use of all_seq_ID (mapping vector) to extract sequences name retained 
   # in full BG that have the same sequence to those inserted by the user. 
   x <- all_seq_ID[ID]
   x <- unique(x) 
   
   # NA removal. These correspond to sequences removed by .clean_sequence()
   rem_names <- names(x[is.na(x)])
   
   if(length(rem_names)> 0){
     warning(paste('Found', length(rem_names), 'sequences with more than 50% of 
                  N content or with a length different from the reference. 
                  Removing the following sequences:', 
                   paste(rem_names, collapse = ", ")))
   }
   
   x <- x[!is.na(x)]
   
   # See ps_scan for details 
   
   pfms <- lapply(
     full_pfms, 
     FUN = ps_scan, 
     x, 
     BG = FALSE,
     use_full_BG = TRUE)
   
   BiocGenerics::do.call(PSMatrixList, pfms)
}

#' Create a Summary Table of PscanR Results  
#' 
#' This function generates a table summarizing the statistical analysis of motif 
#' enrichment, stored in a `PSMatrixList` object. It retrieves pre-computed 
#' background and foreground statistics (the alignment scores between Pscan 
#' input regulatory sequences and position weight matrices) for each PWM. 
#' Then it compares motif occurrences in the foreground (e.g., 
#' coexpressed/coregulated promoter sequences) with the background (e.g., all 
#' promoters in an organism) using Z-scores, p-values, and FDR correction.
#' Returns a table ordered by decreasing `ZSCORE` and increasing `P.VALUE`.  
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` or `pscan_full_bg()` functions. 
#'
#' @return
#' A data.frame with matrices ordered by increasing P.VALUE and decreasing 
#' ZSCORE, prioritizing significant and strongly enriched motifs.
#' It contains the following columns: 
#' \itemize{
#'   \item "NAME": The matrices identifiers.
#'   \item "BG_AVG": The average background score of TFBSs for each PWM.
#'   \item "BG_STDEV": The standard deviation of the background scores of TFBSs
#'   for each PWM.
#'   \item "FG_AVG": The average foreground score (the TFBS scores found by 
#'   scanning the chosen subset of regulatory regions) for each PWM.
#'   \item "ZSCORE": The Z-score for each PWM computed as 
#'   (FG_AVG-BG_AVG)/BG_STDEV. It represent a statistical measure for motif 
#'   enrichment of the scanned sequences in respect to the background.
#'   \item "P.VALUE": The p-value for each PWM.
#'   \item "FDR": The adjusted p-value (False Discovery Rate) calculated using 
#'   the Benjamini-Hochberg correction.
#' }
#' 
#' @details
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
#' 
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' ps_results_table(results)
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

#' Generate a Z-score Table for Motif Hits
#' 
#' This function creates a matrix of Z-scores for motif occurrences across 
#' multiple PSMatrix contained in a `PSMatrixList` object. Each row corresponds 
#' to a sequence, and each column represents a motif. 
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
#' @details
#' The Z-Score represents the statistical significance of the alignment scores 
#' for regulatory sequences relative to background expectation. A high Z-score 
#' suggests strong motif enrichment in the foreground compared to 
#' the background.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @examples
#' # The generation of the example might take few minutes 
#' 
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:10]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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
#'    of `pscan()` or `pscan_full_bg()` functions. 
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to include in the analysis. The default is `0.01`.
#' @param ... Additional user defined arguments to customize the heatmap 
#'    settings, such as color palettes or clustering object.
#'   
#' @details
#' 
#' The heatmap represents the correlation of motif Z-scores, helping to identify 
#' clusters of motifs that show similar enrichment patterns across sequences. 
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
#'    \item `cluster_rows` and `cluster_cols`, set to `TRUE` by default.
#'    \item `color`, which uses a blue-white-red palette.
#'    \item The main `Pscan Score Correlation Heatmap`.
#'    \item TF names are showed as column labels.
#'    }
#'    
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#'
#' @return A heatmap plot showing Z-score correlations for selected 
#'    transcription factors. 
#' 
#' @examples
#' # The generation of the example might take few minutes
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:25]
#'
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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
#' of motif hits based on a specified false discovery rate threshold. 
#' The heatmap helps visualize where significant motif hits occur within the
#' analyzed sequences.
#' 
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated 
#'    metadata (foreground and background statistics). Typically is the output 
#'    of `pscan()` or `pscan_full_bg()` function. 
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to be included in the analysis. The default is set to `0.01`.
#' @param shift Integer. A value to shift the reported positions of motif hits. 
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
#'      with customizable settings.}
#' 
#' Default settings, that can be changed by the users, are:
#' \itemize{
#'   \item `cluster_rows` and `cluster_cols`, set to `TRUE` by default.
#'   \item `color`, which uses a white-yellow-red palette by default.
#'   \item The main `Pscan Hits Position Heatmap`.
#'   }
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @return A heatmap plot showing positional hits distribution for selected 
#'    transcription factors.
#' 
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:25]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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
  
  pos_mat <- matrix(data = NA, nrow = ps_fg_size(pfms[[1]]), 
                    ncol = length(topn))
  
  for(v in topn)
    pos_mat[,v] <- ps_hits_pos(pfms[[row.names(res_table)[v]]], 
                               pos_shift = shift)
  
  colnames(pos_mat) <- res_table$NAME[topn]
  rownames(pos_mat) <- ps_seq_names(pfms[[1]])
  
  do.call(pheatmap, c(list(pos_mat), final_args))
  
  invisible(pos_mat)
}

#' Pscan Density Plot of Hits along Promoters 
#' 
#' This function creates a density plot representing the distribution of hits 
#' along the promoter sequences based on their position and score for a specific 
#' PWM. It helps visualize where motif occurrences are concentrated.
#' 
#' @param pfm A `PSMatrix` object. It is a selected matrix from the 
#'    `PSMatrixList`, result of `pscan` function. 
#' @param shift Integer value specifying the positional shift applied to the hit 
#'    positions. Default is `0`.
#' @param st Score threshold used to filter hits. Can be a numeric value to set 
#'    the threshold directly, or a character:
#'    \itemize{
#'      \item `all`: the threshold is set to `0` (All the hits are evaluated).
#'      \item `loose`: uses the background average score as threshold.
#'      \item `strict`: uses the background average score together with the 
#'      background standard deviation as threshold.}
#'    Default is set to loose.
#'      
#' @return A density plot showing the distribution of hits along the promoter 
#'    region.
#'    
#' @details
#' The function filters motif hits based on a specified threshold and generates 
#' a density plot to show their distribution. The function includes a vertical 
#' dashed line marking the mode (the most frequent position along the promoters).
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#'      
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:25]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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

#' Bubble chart of Pscan Motif Score vs Position
#' 
#' This function creates a Bubble Chart to visualize the relationship
#' between position and score of the identified sites along the promoter 
#' sequences in a `PSMatrix` object. 
#' 
#' @param pfm A `PSMatrix` object. It must be processed by the PscanR 
#'    algorithm.
#' @param bubble_color A character string specifying the color of the bubbles 
#'    (default: `"blue"`).
#'    
#' @return A ggplot2 bubble chart. The x-axis corresponds to the position of 
#' motif hits along the promoter, the y-axis corresponds to the score of motif 
#' hits, and the size of the bubbles represents the count of occurrences for 
#' each score-position combination. 
#' 
#' @details
#' The function computes the count of occurrences for each combination of 
#' PS hits score and position from the input `PSMatrix` and visualizes 
#' this data using a bubble chart. The user can modify the color of bubbles. 
#' Default is `blue`.
#' 
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#' 
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:25]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
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
  
  ggplot(data_sum, aes(x = .data[[colnames(data_sum)[2]]], 
                       y = .data[[colnames(data_sum)[1]]], size = .data$Count)) +
    geom_point(alpha=0.5, color = bubble_color) +
    scale_size_continuous(breaks = sort(unique(data_sum$Count)), 
                          guide = guide_legend(title = "Occurrences")) +
    labs(x = "PS Hits Position", y = "PS Hits Score", 
         title = paste(pfm@name,"Bubble Chart of Score vs Position Hits")) +
    theme_minimal()
  
}

#' Density Plot of Distances between Identified Motif Hits in two PSMatrix 
#' Object
#' 
#' This function visualizes the density plot of distances between identified 
#' hits sites in two `PSMatrix` object. The distance between hits is calculated 
#' for each sequence that is present in both matrices. 
#' It allows to filter the identified sites based on a specified threshold value. 
#' 
#'
#' @param M1 A `PSMatrix` object. It must be processed by the PscanR 
#'    algorithm. 
#' @param M2 A `PSMatrix` object. It must be processed by the PscanR 
#'    algorithm. 
#' @param st1 Score threshold used to filter hits for the first `PSMatrix` . 
#'    Can be a numeric value to set the threshold directly, or a character:
#'    \itemize{
#'      \item `all`: the threshold is set to `0` (All hits are evaluated).
#'      \item `loose`: uses the background average as threshold.
#'      \item `strict`: uses the background average plus the background standard 
#'      deviation as threshold.}
#'    Default is set to `loose`.
#' @param st2 Score threshold used to filter hits for the second `PSMatrix`. 
#'    Can be a numeric value to set the threshold directly, or a character, as 
#'    for st1. 
#'    Default is set to `loose`. 
#'
#' @return A density plot showing the distribution of distances between 
#'    identified motif hits in \code{M1} and \code{M2}. The x-axis represents 
#'    the distances between corresponding hits, and the y-axis represents 
#'    the density of those distances.
#' 
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_bg_std_dev}}, 
#' \code{\link{ps_hits_score}}, \code{\link{ps_hits_pos}}
#' 
#' @details
#' This function uses example datasets located in the `extdata/` directory for 
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the 
#' examples.
#'
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions: 
#' # -200 +50 bp in respect to the TSS. 
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[25:50]
#' 
#' # Retrieve Background PWMs
#' J2020_PSBG <- generate_psmatrixlist_from_background('Jaspar2020', 
#'                                                     'hs', c(-200,50), 'hg38')
#' 
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, J2020_PSBG, 
#'                  BPPARAM = BiocParallel::SnowParam(1))
#' # Use MulticoreParam() for Unix systems (See BiocParallel package).
#' 
#' pfm1 <- results[[1]]
#' pfm2 <- results[[2]]
#' ps_density_distances_plot(pfm1, pfm2, 'all', 'loose')
#' @export
ps_density_distances_plot <- function(M1, M2, st1 = ps_bg_avg(M1), st2 = ps_bg_avg(M2))
{
  if (!is(M1, "PSMatrix") || !is(M2, "PSMatrix"))
    stop("Both object must be of class PSMatrix")
  
  if (is.character(st1)){
    st1 <- switch(st1, "all" = 0, "loose" = ps_bg_avg(M1), 
                  "strict" = ps_bg_avg(M1) + ps_bg_std_dev(M1),
                  {
                    warning("Invalid value for st1, reverting to loose")
                    ps_bg_avg(M1)
                  })
  }
  
  if (is.character(st2)) {
    st2 <- switch(st2, "all" = 0, "loose" = ps_bg_avg(M2),
                  "strict" = ps_bg_avg(M2) + ps_bg_std_dev(M2),
                  {
                    warning("Invalid value for st2, reverting to loose")
                    ps_bg_avg(M2)
                  })
  }
  
  scores1 <- ps_hits_score(M1)
  g_scores1 <- scores1 >= st1
  scores2 <- ps_hits_score(M2)
  g_scores2 <- scores2 >= st2
  
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
