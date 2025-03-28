#' Get the Transcript ID Legend from a PSMatrixList object
#' 
#' This method retrieves the value stored in the `transcriptIDLegend` slot of
#' a `PSMatrixList` object.
#' 
#' @param x An object of class `PSMatrixList`.
#' 
#' @return A named character vector containing transcript ID mappings where 
#'    values are the sequence identifiers maintained during the generation of 
#'    the background by the unique() function, while the names are the 
#'    original sequence identifiers for the organism of study.    
#' 
#' @export
setMethod("transcriptIDLegend", "PSMatrixList", function(x) {
  out <- x@transcriptIDLegend

  return(out)
})

#' Get Background Average Score
#'
#' Retrieves the background average score stored in a `PSMatrix` object.
#' This score represents the average binding score of promoter regions 
#' when scanned with a PWM.
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`, but this argument is currently unused as 
#'    the function returns a single numeric value. It is included for 
#'    consistency with other matrix-related methods.
#' 
#' @return A numeric value representing the background average score. This value 
#'    is computed as the mean of the PWM scores obtained from scanning all 
#'    promoter regions in the organism.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_avg(pfm1) 
#'
#' @export
setMethod("ps_bg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_avg
  
  return(out)
})

#' Get Foreground Average Score
#'
#' Retrieves the foreground average score stored in a `PSMatrix` object. 
#' This score represents the average binding score of a set of promoter regions 
#' of co-expressed or co-regulated genes when scanned with a PWM.
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`, but this argument is currently unused as 
#'    the function returns a single numeric value. It is included for 
#'    consistency with other matrix-related methods.
#' 
#' @return A numeric value representing the foreground average score. This value 
#'    is computed as the mean of the PWM scores obtained from scanning a set of
#'    promoter regions of genes of interest.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_fg_avg(pfm1)
#' 
#' @export
setMethod("ps_fg_avg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_avg
  
  return(out)
})

#' Get Z-Score
#' 
#' Retrieves the Z-Score stored in a `PSMatrix` object.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in 
#'    the output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the `PSMatrix` Z-score. This associates 
#' with each profile the probability of obtaining the same score on a random 
#' sequence set. 
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_zscore(pfm1)
#'
#' @export
setMethod("ps_zscore", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_zscore
  
  return(out)
})

#' Get P-Value
#' 
#' Retrieves the P-Value stored in a `PSMatrix` object. The P-value quantifies
#' the statistical significance of motif enrichment in the input sequences in
#' respect to the background. 
#' 
#' @param x A `PSMatrix` object, typically the result of the `Pscan` algorithm. 
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @return A numeric value representing the statistical significance of motif 
#'     enrichment. Lower values suggest stronger motif enrichment in the 
#'     foreground set.
#'     
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_pvalue(pfm1)
#' 
#' @export
setMethod("ps_pvalue", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_pvalue
  
  return(out)
})

#' Get Matched Oligonucleotide Sequences 
#' 
#' Retrieves the oligonucleotide sequences that match the motif in the 
#' foreground sequences, returning a character vector.
#' 
#' @param x A `PSMatrix` object, typically the result of the `Pscan` algorithm.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A character vector containing the sequences of motif matches 
#'     (oligonucleotides) in the foreground set. 
#'     
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_oligo(pfm1)
#' 
#' 
#' @export
setMethod("ps_hits_oligo", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_oligo

  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get Matched Oligonucleotide Sequences of a Background Dataset
#' 
#' Retrieves the oligonucleotide sequences that match the motif in the 
#' background sequences, returning a character vector.
#' 
#' @param x A `PSMatrix` object, typically the result of the `Pscan` algorithm.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @details
#' This method is specifically designed for background datasets,
#' where motif hit oligonucleotide sequences are precomputed for all 
#' promoter sequences. 
#' The function extracts the character values from the ps_hits_oligo_bg slot 
#' and, if applicable, assigns sequence names using .ps_bg_seq_names().
#' 
#' These background metrics are particularly useful in motif enrichment analyses, 
#' as they allow pscan() to compare foreground promoter sequences against a 
#' reference distribution without the need for recomputation.
#'    
#' @return A character vector containing the sequences of motif matches 
#'     (oligonucleotides) in the banckground set. 
#'     
#' @examples
#' full_pfm1_path <- system.file("extdata", "full_pfm1.rds", package = "PscanR")
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_oligo_bg(full_pfm1)
#' 
#' @export
setMethod("ps_hits_oligo_bg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_oligo_bg
  
  out <- .ps_bg_seq_names(x, out)
  
  return(out)
})

setMethod(".ps_bg_seq_names", "PSMatrix", function(x, out) {
  
  if(!any(is.na(x@ps_bg_seq_names)))
    names(out) <- x@ps_bg_seq_names
  
  return(out)
})

setMethod(".ps_seq_names", "PSMatrix", function(x, out) {
  
  if(!any(is.na(x@ps_seq_names)))
    names(out) <- x@ps_seq_names
  
  return(out)
})

#' Get Background Standard Deviation
#' 
#' Retrieves the background standard deviation of PWM scores stored in a 
#' `PSMatrix` object. This value quantifies the variability in the background 
#' motif binding scores.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A numeric value representing the standard deviation of the background 
#'    PWM scores.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_std_dev(pfm1)
#' 
#' @export
setMethod("ps_bg_std_dev", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_std_dev
  
  return(out)
})

#' Get Background Size
#' 
#' Retrieves the number of promoter regions used to compute the background 
#' statistics in a `PSMatrix` object. This represents the total number of 
#' promoters expressed in the organism of study. 
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer value representing the total number of promoter regions
#'    used for background scoring. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_size(pfm1)
#' 
#' @export
setMethod("ps_bg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_size
  
  return(out)
})

#' Get Foreground Size
#' 
#' Retrieves the number of promoter sequences given as input to the Pscan 
#' algorithm in a `PSMatrix` object. This represents the total number of 
#' sequences analyzed for motif enrichment.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer representing the total number of promoter sequences used 
#'    as input to Pscan.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_fg_size(pfm1)
#' 
#' @export
setMethod("ps_fg_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_fg_size
  
  return(out)
})

#' Compute Hits Size
#' 
#' Retrieves the number of motif hits detected by the Pscan algorithm in a 
#' `PSMatrix` object. A motif hit represents a significant match between a 
#' promoter sequence and the Position Weight Matrix (PWM).
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return An integer representing the total number of motif hits detected in 
#'    the input promoter sequences.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_size(pfm1)
#' 
#' @export
setMethod("ps_hits_size", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- length(x@ps_hits_pos)
  
  return(out)
})

#' Get Hits Score
#' 
#' Retrieves the motif hit scores for each promoter sequence in a `PSMatrix` 
#' object.
#' These scores represent the binding affinity or enrichment level of promoter 
#' sequences when scanned with a Position Weight Matrix (PWM).
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A named numeric vector where names correspond to promoter sequence 
#'    identifiers and values represent their respective motif hit scores.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_score(pfm1)
#' 
#' @export
setMethod("ps_hits_score", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_score
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get Hits Score for Background Dataset
#' 
#' Retrieves the motif hit scores for each promoter sequence in a `PSMatrix` 
#' object computed on the background dataset (all the promoter in a specific
#' organism).
#' These scores represent the binding affinity or enrichment level of promoter 
#' sequences when scanned with a Position Weight Matrix (PWM).
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @details
#' This method is specifically designed for background datasets,
#' where motif hit scores are precomputed for all promoter sequences. 
#' The function extracts the scores from the ps_hits_score_bg slot and, 
#' if applicable, assigns sequence names using .ps_bg_seq_names().
#' 
#' These background scores are particularly useful in motif enrichment analyses, 
#' as they allow pscan() to compare foreground promoter sequences against a 
#' reference distribution without the need for recomputation.  
#'
#' @return A named numeric vector where names correspond to promoter sequence 
#'    identifiers and values represent their respective motif hit scores.
#' 
#' @examples
#' full_pfm1_path <- system.file("extdata", "full_pfm1.rds", package = "PscanR")
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_score_bg(full_pfm1)
#' 
#' @export
setMethod("ps_hits_score_bg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_score_bg
  
  out <- .ps_bg_seq_names(x, out)
  
  return(out)
})

#' Get Motif Hit Z-score
#' 
#' Computes the Z-Scores for motif hit scores in a `PSMatrix` object. 
#' The Z-score measures how unusual a motif score is compared to background 
#' sequences (all the promoters expressed in an organism). Higher Z-scores mean 
#' stronger motif enrichment, suggesting potential regulatory significance.
#'
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A named numeric vector where names correspond to promoter sequence 
#'    identifiers and values represent their respective Z-scores.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_z(pfm1)
#' 
#' @export
setMethod("ps_hits_z", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- (x@ps_hits_score - x@ps_bg_avg) / x@ps_bg_std_dev
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get Motif Hit Strand Information
#' 
#' Retrieves the strand information (`+` or `-`) of motif hits in a 
#' `PSMatrix` object. This indicates whether a motif was detected on the forward 
#' (`+`) or reverse (`-`) strand of the promoter sequence.
#' The Pscan algorithm scans both the strands of the promoter sequences to 
#' ensure that no potential binding sites are missed.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A named character vector where names correspond to promoter sequence 
#'    identifiers, and values represent the strand (`+` or `-`) on which the 
#'    motif was detected.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_strand(pfm1)
#' 
#' @export
setMethod("ps_hits_strand", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_strand
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})


#' Get Motif Hit Strand Information of a Background Dataset
#' 
#' Retrieves the strand information (`+` or `-`) of motif hits in a 
#' `PSMatrix` object computed on the background dataset (all the promoter in a 
#' specific organism). This indicates whether a motif was detected on 
#' the forward (`+`) or reverse (`-`) strand of the promoter sequence.
#' The Pscan algorithm scans both the strands of the promoter sequences to 
#' ensure that no potential binding sites are missed.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @details
#' This method is specifically designed for background datasets,
#' where motif hit strands are precomputed for all promoter sequences. 
#' The function extracts the strand values from the ps_hits_strand_bg slot and, 
#' if applicable, assigns sequence names using .ps_bg_seq_names().
#' 
#' These background metrics are particularly useful in motif enrichment analyses, 
#' as they allow pscan() to compare foreground promoter sequences against a 
#' reference distribution without the need for recomputation.
#'    
#' @return A named character vector where names correspond to promoter sequence 
#'    identifiers, and values represent the strand (`+` or `-`) on which the 
#'    motif was detected.
#' 
#' @examples
#' full_pfm1_path <- system.file("extdata", "full_pfm1.rds", package = "PscanR")
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_strand_bg(full_pfm1)
#' 
#' @export
setMethod("ps_hits_strand_bg", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_hits_strand_bg
  
  out <- .ps_bg_seq_names(x, out)
  
  return(out)
})

#' Get Motif Hit Positions
#' 
#' Retrieves the positions of hits stored in a `PSMatrix` object. These 
#' positions indicate where the motifs are located in each promoter sequence. 
#' The positions can be optionally shifted by a specified value.
#' 
#' @param x A `PSMatrix` object.
#' @param pos_shift Integer. Specifies the amount to shift the position. 
#'    Default is set to `0`.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A named integer vector where names correspond to promoter sequence 
#'    identifiers, and values represent the adjusted motif hit positions.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_pos(pfm1)
#' 
#' @export
setMethod("ps_hits_pos", "PSMatrix", function(x, pos_shift = 0L, 
                                              withDimnames = TRUE) {
  out <- x@ps_hits_pos + pos_shift
  
  out <- .ps_seq_names(x, out)
  
  return(out)
})

#' Get Motif Hit Positions in a Backrgound Dataset
#' 
#' Retrieves the positions of hits stored in a `PSMatrix` object computed on 
#' the background dataset (all the promoter in a specific organism). These 
#' positions indicate where the motifs are located in each promoter sequence. 
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @details
#' This method is specifically designed for background datasets,
#' where motif hit positions are precomputed for all promoter sequences. 
#' The function extracts the position values from the ps_hits_pos_bg slot and, 
#' if applicable, assigns sequence names using .ps_bg_seq_names().
#' 
#' These background metrics are particularly useful in motif enrichment analyses, 
#' as they allow pscan() to compare foreground promoter sequences against a 
#' reference distribution without the need for recomputation.
#' 
#' @return A named integer vector where names correspond to promoter sequence 
#'    identifiers, and values represent the motif hit positions.
#' 
#' @examples
#' full_pfm1_path <- system.file("extdata", "full_pfm1.rds", package = "PscanR")
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_pos_bg(full_pfm1)
#' 
#' @export
setMethod("ps_hits_pos_bg", "PSMatrix", function(x,withDimnames = TRUE) {
  out <- x@ps_hits_pos_bg 
  
  out <- .ps_bg_seq_names(x, out)
  
  return(out)
})

#' Get the Sequences Name 
#' 
#' Retrieves the names or identifiers of the promoter sequences in a `PSMatrix` 
#' object. These names correspond to the promoter regions that were scanned by 
#' the Pscan algorithm.
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'    
#' @return A character vector of names.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_seq_names(pfm1)
#' 
#' @export
setMethod("ps_seq_names", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_seq_names
  
  return(out)
})

#' Get the Sequences Identifiers of the Background Dataset 
#' 
#' Retrieves the names or identifiers of the promoter sequences in a `PSMatrix` 
#' object computed on the background dataset (all the promoter in a specific
#' organism).
#' 
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the 
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#' 
#' @details
#' This method is specifically designed for background datasets.
#' The function extracts the sequence identifiers from the ps_bg_seq_names slot 
#' and, if applicable, assigns sequence names using .ps_bg_seq_names().
#' 
#' @return A character vector of names.
#' 
#' @examples
#' full_pfm1_path <- system.file("extdata", "full_pfm1.rds", package = "PscanR")
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_bg_seq_names(full_pfm1)
#' 
#' @export
setMethod("ps_bg_seq_names", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@ps_bg_seq_names
  
  return(out)
})


setMethod(".PS_PSEUDOCOUNT", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_PSEUDOCOUNT
  
  return(out)
})


setMethod(".PS_ALPHABET", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@.PS_ALPHABET
  
  return(out)
})

setMethod("all_sequences_ID", "PSMatrix", function(x, withDimnames = TRUE) {
  out <- x@all_sequences_ID
  
  return(out)
})

#' Get Table of Motif Hits
#' 
#' Creates a data frame summarizing the motif hits from a `PSMatrix` object, 
#' including the motif hit score, position, strand, and corresponding oligo 
#' sequence.
#' The resulting table is sorted by motif score (descending). When the score is 
#' the same between two sequences, they are sorted by ascending position.
#' 
#' @param x An object of class `PSMatrix`. Should contain the following 
#'    slots:
#'    \itemize{
#'      \item `ps_hits_score`: a numeric vector of hit scores.
#'      \item `ps_hits_strand`: a character vector of strand information
#'      (`-` and `+`).
#'      \item `ps_hits_oligo`: a character vector of oligo sequences.}
#' @param pos_shift Integer. Value for which the position gets shifted.
#'    Default is `0`, meaning no shift.
#' @param withDimnames Logical, whether to include dimension names in the output, 
#'    if they exist in the object.
#'    Default set to `TRUE`.
#'  
#' @return A `data.frame` ordered by decreasing score value, 
#' with the following columns:
#' \itemize{
#'   \item `SCORE`: the motif hit score
#'   \item `POS`: the position of the motif hit
#'   \item `STRAND`: the strand orientation for each hit
#'   \item `OLIGO`: the oligo sequence corresponding to each hit} 
#' Row names correspond to the sequences name. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_table(pfm1)
#' 
#' @export
setMethod("ps_hits_table", "PSMatrix", function(x, pos_shift = 0L, 
                                                withDimnames = TRUE) {
  
  out <- data.frame("SCORE" = x@ps_hits_score, 
                    "POS" = ps_hits_pos(x, pos_shift = pos_shift), 
                    "STRAND" = x@ps_hits_strand,
                    "OLIGO" = DNAStringSet(x@ps_hits_oligo), 
                    row.names = x@ps_seq_names)
  
  out <- out[with(out, order(SCORE, POS, decreasing = c(TRUE,FALSE))),]
  
  return(out)
})


setMethod(".ps_add_hits", "PSMatrix", 
          function(x, Pos, Strand, Score, Oligo, BG = FALSE, use_full_BG = FALSE,
                   fullBG = FALSE, withDimnames = TRUE) {
  
  x@ps_hits_pos <- Pos
  x@ps_hits_strand <- Strand
  x@ps_hits_score <- Score
  if (!use_full_BG)
    x@ps_hits_score <- .ps_norm_score(x)
  
  if(BG)
  {
    ps_bg_size(x) <- length(x@ps_hits_pos)
    ps_bg_avg(x) <- mean(x@ps_hits_score, na.rm = TRUE)
    ps_bg_std_dev(x) <- sd(x@ps_hits_score, na.rm = TRUE)
    if(ps_bg_std_dev(x) == 0)
      ps_bg_std_dev(x) <- 0.00001
    
    if(fullBG){
      x@ps_hits_pos_bg <- Pos
      x@ps_hits_strand_bg <- Strand
      x@ps_hits_score_bg <- x@ps_hits_score
      names(x@ps_hits_score_bg) <- x@ps_bg_seq_names
      x@ps_hits_oligo_bg <- Oligo
    }
    
    x@ps_hits_pos <- integer()
    x@ps_hits_strand <- character()
    x@ps_hits_score <- numeric()
  }
  else {
    if(!is.na(x@ps_bg_avg) && !is.na(x@ps_bg_std_dev))
    {
      ztest <- z.test(x@ps_hits_score, 
                      mu = x@ps_bg_avg, 
                      sigma.x = x@ps_bg_std_dev, 
                      alternative = "greater")
      
      x@ps_zscore <- ztest$statistic["z"]
      x@ps_pvalue <- as.numeric(ztest$p.value)
      x@ps_fg_avg <- mean(x@ps_hits_score, na.rm = TRUE)
      x@ps_fg_size <- length(x@ps_hits_pos)
      x@ps_hits_oligo <- Oligo
    }
    
    x@ps_hits_pos_bg <- integer()
    x@ps_hits_strand_bg <- character()
    x@ps_hits_score_bg <- numeric()
    x@ps_hits_oligo_bg <- character()
  }
  
  return(x)
})


setMethod(".ps_norm_score", "PSMatrix", function(x) {
  
  ps_score <- 1 + ((maxScore(Matrix(x)) - ps_hits_score(x)) / 
         (minScore(Matrix(x)) - maxScore(Matrix(x))))
  
  return(ps_score)
})


setMethod(".ps_bg_from_table", "PSMatrix", function(x, short.matrix) {
  
  if(any(row.names(short.matrix) == ID(x)))
  {
    x@ps_bg_size <- as.integer(short.matrix[ID(x),"BG_SIZE"])
    x@ps_bg_avg <-  as.numeric(short.matrix[ID(x),"BG_MEAN"])
    x@ps_bg_std_dev <-  as.numeric(short.matrix[ID(x),"BG_STDEV"])
  }
  else
  {
    warning(paste("No background values found for", ID(x), name(x)))
  }

  return(x)
})


#' @importMethodsFrom TFBSTools Matrix
setMethod(".ps_norm_matrix", "PSMatrix", function(x){
  
  mx <- Matrix(x)
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
  mx <- mx + x@.PS_PSEUDOCOUNT
  
  mx <- sweep(mx, 2, colSums(mx), FUN = "/")
  
  mx <- log(mx)
  
  Matrix(x) <- mx
  
  return(x)
  
})

#' Perform a Scan of DNA Sequences Using a `PSMatrix` Object
#' 
#' This method performs a scan of DNA sequences using a `PSMatrix` object 
#' to identify potential hits based on the matrix score. The scan is 
#' performed on both strands (the direct and reverse complement strands) of the 
#' sequences to ensure that no potential binding sites are missed. 
#' 
#' @param x A `PSMatrix` object representing the motif and related statistics 
#'     to be used for scanning DNA sequences.
#' @param seqs `A DNAstringSet` object representing the sequences to be
#'    scanned for motif occurrences.
#' @param BG A logical value indicating whether to calculate background 
#'    statistics.
#'    Default is set to `FALSE`.
#' @param use_full_BG A logical value (default is `FALSE`). If `TRUE`, the method 
#'    assumes that the PSMatrix represent a special "full-background" set. 
#'    In this case, the `seqs` parameter should be a vector of sequence names 
#'    (instead of actual DNA sequences). The function handles this by matching 
#'    sequence names against precomputed background data in the `PSMatrix`.
#' @param fullBG Logical, default is `FALSE`. If `TRUE` it computes a series 
#'    of computation for the generation of a "full-background", a special case
#'    of background that retains all the background hits score, position, 
#'    strand, and oligonucleotide sequence for each regulatory sequence scanned 
#'    with a PWM.  
#' 
#' @return A `PSMatrix` object, with updated information about the motif hits 
#'     in the sequences. This includes the positions, strands, scores, and 
#'     oligos (sequence motifs) where the hits occurred.
#' 
#' @details
#' 
#' This method is designed to handle different types of sequence scanning 
#' scenarios. It can handle both regular DNA sequences (using a `DNAStringSet`) 
#' and specialized background sequence sets (either using a background flag 
#' or a full-background flag). 
#' 
#' If `use_full_BG` is set to `TRUE`, the function assumes that `seqs` contains 
#' sequence identifiers rather than the sequences themselves. In this case, the 
#' method matches the sequence names to those in the `PSMatrix`'s background hit 
#' data and retrieves the corresponding binding information (score, strand, 
#' position, oligo).
#' 
#' When `use_full_BG` is set to `FALSE`, the function scan both the forward and 
#' reverse complement strands of the sequences to ensure all potential binding 
#' sites are detected. Optionally, the background statistics (background average 
#' and standard deviation) can be computed and used during scanning when `BG` is 
#' set to `TRUE`.
#' 
#' @examples
#' BG_matrices <- generate_psmatrixlist_from_background('Jaspar2020', 'hs', 
#'                                                      c(-200,50),'hg38')
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)
#' prom_seq <- prom_seq[1:25]
#' 
#' scanned_result <- ps_scan(BG_matrices[[3]], prom_seq)
#' scanned_result
#' 
#' @export
setMethod("ps_scan", "PSMatrix", function(x, seqs, BG = FALSE, 
                                          use_full_BG = FALSE, fullBG = FALSE){
  
  if(!is(seqs, "DNAStringSet") && !use_full_BG)
    stop("seqs is not an object of DNAStringSet class")
  
  if (BG || use_full_BG) {
    if(use_full_BG)
      x@ps_bg_seq_names <- seqs
    else 
      x@ps_bg_seq_names <- names(seqs)
  } else {
    x@ps_seq_names <- names(seqs)
  }
  
  seqs <- as.character(seqs)
  
  if(use_full_BG == TRUE){
    indices <- match(seqs, sub("\\..*$", "", names(x@ps_hits_score_bg)))
    
    res <- list(
      score = x@ps_hits_score_bg[indices],
      strand = x@ps_hits_strand_bg[indices],
      pos = x@ps_hits_pos_bg[indices],
      oligo = x@ps_hits_oligo_bg[indices]
    )
    
    x <- .ps_add_hits(x, Score = res$score, Strand = res$strand, Pos = res$pos, 
                      Oligo = res$oligo, BG = BG, use_full_BG = use_full_BG)
    x@ps_bg_seq_names <- character()
    
  } else {
    
    rc_x <- reverseComplement(x)
    
    Margs <- list(numx = as.numeric(Matrix(x)), 
                numx_rc = as.numeric(Matrix(rc_x)),
                ncolx = (0:(ncol(Matrix(x)) - 1))*length(.PS_ALPHABET(x)), 
                AB = .PS_ALPHABET(x)) 
    
    res <- mapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs)
    
    x <- .ps_add_hits(x, Score = as.numeric(res["score",]), 
                      Strand = as.character(res["strand",]), 
                      Pos = as.integer(res["pos",]), 
                      Oligo = as.character(res["oligo",]), BG = BG, 
                      use_full_BG = use_full_BG, fullBG = fullBG)
    if(!fullBG)
      x@ps_bg_seq_names <- character()
  }
  
  return(x)
  
})


#' @importMethodsFrom Biostrings maxScore minScore
setMethod(".ps_scan_s", "PSMatrix", function(x, Seq, numx, numx_rc, ncolx, AB){
  
  subS <- strsplit(substring(Seq, seq_len((nchar(Seq) - length(x) + 1)), 
                             length(x):nchar(Seq)),"",
                   fixed = TRUE)
  prot <- numeric(1)
  
  scores <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, 
                   x = numx, AB = AB, ncolx = ncolx)
  scores_rc <- vapply(subS, FUN = .ps_assign_score, FUN.VALUE = prot, 
                      x = numx_rc, AB = AB, ncolx = ncolx)
  
  mscore_pos <- which.max(scores)
  mscore_rc_pos <- which.max(scores_rc)
  
  res <- list(score = numeric(), strand = character(), pos = integer(), 
              oligo = character())
  
  if(scores[mscore_pos] >= scores_rc[mscore_rc_pos])
  {
    res$score <- scores[mscore_pos]
    res$strand <- "+"
    res$pos <- mscore_pos
    res$oligo <- paste(subS[[mscore_pos]], collapse = '')
  }
  else
  {
    res$score <- scores_rc[mscore_rc_pos]
    res$strand <- "-"
    res$pos <- mscore_rc_pos
    res$oligo <- paste(subS[[mscore_rc_pos]], collapse = '')
  }
  
  return(res)
})


#' @importMethodsFrom Biostrings reverseComplement 

#setMethod(".ps_assign_score", "PSMatrix", function(x, S){
#  sum(Matrix(x)[matrix(data = c(.PS_ALPHABET(x)[S], 1:length(x)), ncol = 2, nrow = length(x))])
#})

.ps_assign_score <- function(S, x, AB,ncolx)
{
  sum(x[ncolx+AB[S]]) #Assign score to oligo
}

#' Validate a PSMatrix object
#' 
#' This function checks if the input `PSMatrix` object properties meet the ones
#' desired.
#' 
#' @param object A `PSMatrix` object.
#' 
#' @details
#' The function `validPSMatrix` ensures that the `PSMatrix` object has correctly 
#' formatted background and foreground averages, standard deviations, 
#' and hit-related values. 
#' It performs the following steps:
#' \itemize{
#'   \item `ps_bg_avg`,`ps_fg_avg`, and `ps_bg_std_dev` must be of length 1.
#'   \item `ps_bg_avg` and `ps_bg_std_dev` value must be between 0 and 1, 
#'   excluding 0 for `ps_bg_std_dev`.
#'   \item the length of `ps_hits_pos`, `ps_hits_strand`, and `ps_hits_score` 
#'   vectors must be equal.}
#' 
#' @return If all the checks are satisfied, returns `TRUE`. Otherwise, a
#' string describing the reason of failure. 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' validPSMatrix(pfm1)
#' 
#' @export
validPSMatrix <- function(object)
{
  if(length(object@ps_bg_avg) != 1)
    return("Background average must be of length 1")
  if(length(object@ps_fg_avg) != 1)
    return("Foreground average must be of length 1")
  if(length(object@ps_bg_std_dev) != 1)
    return("Background stdev must be of length 1")
  if((object@ps_bg_avg < 0 || object@ps_bg_avg > 1) && !is.na(object@ps_bg_avg)) 
    return(paste("Invalid value for Background average: ", object@ps_bg_avg))
  if((object@ps_bg_std_dev < 0 || object@ps_bg_std_dev > 1) && !is.na(object@ps_bg_std_dev))
    return(paste("Invalid value for Background stddev: ", object@ps_bg_std_dev))
 # if(object@ps_bg_size < 1000 && !is.na(object@ps_bg_size))
  #  return(paste("Invalid value for Background size: ", object@ps_bg_size, " Background must be of at least 1000 sequences"))
  if(length(object@ps_hits_pos) != length(object@ps_hits_strand) || length(object@ps_hits_pos) != length(object@ps_hits_score))
    return(paste("Invalid PSMatrix object: different values for hits, strands and scores vectors"))
  
  TRUE
}


setValidity("PSMatrix", validPSMatrix)

#' Display Details of a `PSMatrix` object
#' 
#' This method displays a summary of a `PSMatrix` object, including 
#' background and foreground statistics, standard deviation, background 
#' and foreground size, z-scores and P-Value. 
#' 
#' @param object An object of class `PSMatrix`.
#' 
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_fg_avg}}, 
#'    \code{\link{ps_bg_std_dev}}, \code{\link{ps_bg_size}}, 
#'    \code{\link{ps_fg_size}}, \code{\link{ps_zscore}}, \code{\link{ps_pvalue}}
#'    
#' @return Prints a summary of the `PSMatrix` object, including:
#' \itemize{
#'   \item Pscan Background Average
#'   \item Pscan Foreground (the sample) Average
#'   \item Pscan Background Standard Deviation
#'   \item Pscan Background Size
#'   \item Pscan Foreground (the sample) Size
#'   \item Pscan Z-score
#'   \item Pscan p-value (Z-test)
#' }
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' show(pfm1)
#' 
#' @export
#' @importMethodsFrom methods show
setMethod("show", "PSMatrix", function(object) {
  
  callNextMethod()
  
  cat(
      "\nPscan Background Average: ", ps_bg_avg(object),
      "\nPscan Foreground (your sample) Average: ", ps_fg_avg(object),
      "\nPscan Backgroun Stdev: ", ps_bg_std_dev(object),
      "\nPscan Background Size: ", ps_bg_size(object),
      "\nPscan Foreground (your sample) Size: ", ps_fg_size(object),
      "\nPscan Zscore: ", ps_zscore(object),
      "\nPscan pvalue (z-test): ", ps_pvalue(object),
      sep = ""
  )
})

#' Set Background Average 
#' 
#' This method specifically set the background average for an object of 
#' class `PSMatrix`.
#' 
#' @param x An object of class `PSMatrix`.
#' @param value Numeric. The mean background value to be set.
#' 
#' @return The modified `PSMatrix` object with the updated background average. 
#' 
#' @seealso \code{\link{ps_bg_avg}} to retrieve the background average value,
#'    and \code{\link{validObject}} to see which checks are performed on the
#'    input `PSMatrix` object.
#' 
#' @details
#' The background size is the total number of promoters used to compute 
#' background statistics, such as the average score and standard deviation. 
#' This method allows user to modify this value. 
#' 
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_avg<-`(pfm1, value = 0.73473) 
#' 
#' @export
setReplaceMethod("ps_bg_avg", "PSMatrix", function(x,value){
  
  x@ps_bg_avg <- value
  validObject(x)
  x
})

#' Set Background Standard Deviation 
#' 
#' This method updated the background standard deviation stored in a `PSMatrix`
#' object.
#' 
#' @param x A `PSMatrix` object.
#' @param value Numeric. The value representing the new background 
#'    standard deviation.
#' 
#' @return The modified `PSMatrix` object with the updated background 
#'    standard deviation. 
#' 
#' @seealso \code{\link{ps_bg_std_dev}} to retrieve the background 
#'    standard deviation, and \code{\link{validObject}} to see which checks 
#'    are performed on the input `PSMatrix` object.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_std_dev<-`(pfm1, value = 0.08568925)
#' 
#' @export
setReplaceMethod("ps_bg_std_dev", "PSMatrix", function(x,value){
  
  x@ps_bg_std_dev <- value
  validObject(x)
  x
})

#' Set Background Size in a PSMatrix
#' 
#' This method specifically set the background size for a `PSMatrix` object. 
#' The background size represents the number of promoter regions used in the 
#' background model.
#' 
#' @param x A `PSMatrix` object. 
#' @param value An integer representing the new background size. 
#' 
#' @return The modified `PSMatrix` object with the updated background size. 
#' 
#' @seealso \code{\link{ps_bg_size}} to retrieve the background size, 
#'    and \code{\link{validObject}} to see which checks are performed on the
#'    input `PSMatrix` object.
#' 
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_size<-`(pfm1, value = 25629L)
#' 
#' @export
setReplaceMethod("ps_bg_size", "PSMatrix", function(x,value){
  
  x@ps_bg_size <- value
  validObject(x)
  x
})

#' Convert PFMatrix to PSMatrix
#' 
#' @return A `PSMatrix` object created from the input `PFMatrix` object.
#' 
#' @name PSMatrix
#' @aliases setAs,PFMatrix,PSMatrix-method
#' 
#' @export
#' @importFrom TFBSTools PFMatrix
setAs("PFMatrix", "PSMatrix", function(from){
  
  #.ps_norm_matrix(new("PSMatrix", from, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
  #                    ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01))
  
  PSMatrix(from)
  
})

#PSMatrix <- function(pfm, ps_bg_avg = as.numeric(NA), ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA), 
#                     ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01, ...)

#' Convert PFMatrixList to PSMatrixList
#' 
#' This method converts a `PFMatrixList` object into a `PSMatrixList` object by 
#' applying the conversion for each `PFMatrix` in the list to a `PSMatrix`. 
#' The resulting `PSMatrixList` is created by combining the converted `PSMatrix` 
#' objects into a list and returning the new object.
#' 
#' @return A `PSMatrixList` object, which is a list containing the `PSMatrix` 
#'    objects converted from the `PFMatrixList`.
#' 
#' @name PSMatrixList
#' @aliases setAs,PFMatrixList,PSMatrixList-method
#' 
#' @export
#' @importFrom TFBSTools PFMatrixList
setAs("PFMatrixList", "PSMatrixList", function(from){
  
  to <- lapply(from, as, "PSMatrix")
  
  do.call(PSMatrixList, to)
  
})

