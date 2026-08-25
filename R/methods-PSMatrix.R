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
#'
#' @examples
#' full_pfms_path <- system.file("extdata",
#'     "full_pfms.rds",
#'     package = "PscanR"
#' )
#' full_pfms <- readRDS(full_pfms_path)
#' head(transcriptIDLegend(full_pfms))
#'
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
#' @return A numeric value representing the background average score.
#'    This value is computed as the mean of the PWM scores obtained from
#'    scanning all promoter regions in the organism.
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
#' @return A numeric value representing the foreground average score.
#'    This value is computed as the mean of the PWM scores obtained from
#'    scanning a set of promoter regions of genes of interest.
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
#' @return A numeric value representing the `PSMatrix` Z-score (z-statistic)
#' computed during motif scanning.
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
#' @section Strand:
#' The sequence returned is always the **forward strand** of the promoter,
#' whatever strand the motif matched on, following the original Pscan
#' implementation. Use \code{\link{ps_hits_strand}} to find which hits matched
#' on \code{"-"} and reverse-complement those before comparing them with the
#' motif consensus.
#'
#' @seealso \code{\link{ps_hits_strand}}, \code{\link{ps_hits_pos}},
#'   \code{\link{ps_hits_table}}
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_oligo(pfm1)
#'
#' # A reverse-strand hit reads as the motif only after flipping it.
#' table(ps_hits_strand(pfm1))
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
#' These background metrics are particularly useful in motif enrichment
#' analyses, as they allow \code{pscan()} to compare promoter
#' sequences against a reference distribution without the need for
#' recomputation.
#'
#' @return A character vector containing the sequences of motif matches
#'     (oligonucleotides) in the background set.
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'     "full_pfm1.rds",
#'     package = "PscanR"
#' )
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
    if (!any(is.na(x@ps_bg_seq_names))) {
        names(out) <- x@ps_bg_seq_names
    }

    return(out)
})

setMethod(".ps_seq_names", "PSMatrix", function(x, out) {
    if (!any(is.na(x@ps_seq_names))) {
        names(out) <- x@ps_seq_names
    }

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
#' @return A numeric value representing the standard deviation of the
#'    background PWM scores.
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
#' Retrieves the background promoter region size used to compute the
#' background statistics in a `PSMatrix` object (e.g., 250L).
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'
#' @return An integer value representing the background promoter region size
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
#' @param withDimnames Logical, whether to include dimension names in the
#'    output, if they exist in the object.
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
#' object computed on the background dataset (all promoters in a specific
#' organism).
#' These scores represent the binding affinity or enrichment level of promoter
#' sequences when scanned with a Position Weight Matrix (PWM).
#'
#' @param x A `PSMatrix` object.
#' @param withDimnames Logical, whether to include dimension names in the
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'
#' @details
#' This method is specifically designed for background datasets,
#' where motif hit scores are precomputed for all promoter sequences.
#' The function extracts the scores from the \code{ps_hits_score_bg} slot and,
#' if applicable, assigns sequence names using \code{.ps_bg_seq_names()}.
#'
#' These background scores are particularly useful in motif enrichment
#' analyses, as they allow \code{pscan()} to compare foreground promoter
#' sequences against a reference distribution without the need for
#' recomputation.
#'
#' @return A named numeric vector where names correspond to promoter sequence
#'    identifiers and values represent their respective motif hit scores.
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'     "full_pfm1.rds",
#'     package = "PscanR"
#' )
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
#' sequences (all promoters expressed in an organism). Higher Z-scores mean
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
#' `PSMatrix` object. This indicates whether a motif was detected on the
#' forward (`+`) or reverse (`-`) strand of the promoter sequence.
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
#' `PSMatrix` object computed on the background dataset (all promoters in a
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
#' The function extracts the strand values from the \code{ps_hits_strand_bg}
#' slot and, if applicable, assigns sequence names using
#' \code{.ps_bg_seq_names()}.
#'
#' These background metrics are particularly useful in motif enrichment
#' analyses, as they allow pscan() to compare foreground promoter sequences
#' against a reference distribution without the need for recomputation.
#'
#' @return A named character vector where names correspond to promoter sequence
#'    identifiers, and values represent the strand (`+` or `-`) on which the
#'    motif was detected.
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'     "full_pfm1.rds",
#'     package = "PscanR"
#' )
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
#' The positions can be shifted by a specified value to find the corresponding
#' position in respect to the TSS.
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
#' ps_hits_pos(pfm1, pos_shift = -200) # in respect to the TSS
#'
#' @export
setMethod(
    "ps_hits_pos",
    "PSMatrix",
    function(x, pos_shift = 0L, withDimnames = TRUE) {
    out <- x@ps_hits_pos + pos_shift

    out <- .ps_seq_names(x, out)

    return(out)
})

#' Get Motif Hit Positions in a Background Dataset
#'
#' Retrieves the positions of hits stored in a `PSMatrix` object computed on
#' the background dataset (all promoters in a specific organism). These
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
#' These background metrics are particularly useful in motif enrichment
#' analyses, as they allow pscan() to compare foreground promoter sequences
#' against a reference distribution without the need for recomputation.
#'
#' @return A named integer vector where names correspond to promoter sequence
#'    identifiers, and values represent the motif hit positions.
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'     "full_pfm1.rds",
#'     package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_pos_bg(full_pfm1)
#'
#' @export
setMethod("ps_hits_pos_bg", "PSMatrix", function(x, withDimnames = TRUE) {
    out <- x@ps_hits_pos_bg

    out <- .ps_bg_seq_names(x, out)

    return(out)
})

#' Get Sequence Names
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

#' Get Sequence Identifiers of the Background Dataset
#'
#' Retrieves the names or identifiers of the promoter sequences in a `PSMatrix`
#' object computed on the background dataset (all promoters in a specific
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
#' full_pfm1_path <- system.file("extdata",
#'     "full_pfm1.rds",
#'     package = "PscanR"
#' )
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
#' @param pos_shift Integer. Value for which the position gets shifted in
#' respect
#'    to the TSS.
#'    Default is `0`, meaning no shift.
#' @param withDimnames Logical, whether to include dimension names in the
#'    output, if they exist in the object.
#'    Default set to `TRUE`.
#'
#' @return A `data.frame` ordered by decreasing score value,
#' with the following columns:
#' \itemize{
#'   \item `SCORE`: the motif hit score
#'   \item `POS`: the position of the motif hit
#'   \item `STRAND`: the strand orientation for each hit
#'   \item `OLIGO`: the oligo sequence corresponding to each hit (a
#'   `DNAStringSet`).}
#' Row names correspond to the sequence names.
#'
#' @section Strand:
#' `OLIGO` is the forward strand of the promoter whatever `STRAND` says, as in
#' the original Pscan. A row with `STRAND == "-"` matches the motif only after
#' its oligo is reverse-complemented. See \code{\link{ps_hits_oligo}}.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_table(pfm1)
#'
#' @export
setMethod("ps_hits_table", "PSMatrix", function(x, pos_shift = 0L,
                                                withDimnames = TRUE) {
    out <- data.frame(
        "SCORE" = x@ps_hits_score,
        "POS" = ps_hits_pos(x, pos_shift = pos_shift),
        "STRAND" = x@ps_hits_strand,
        "OLIGO" = DNAStringSet(x@ps_hits_oligo),
        row.names = x@ps_seq_names
    )

    out <- out[with(out, order(SCORE, POS, decreasing = c(TRUE, FALSE))), ]

    return(out)
})


setMethod(
    ".ps_add_hits", "PSMatrix",
    function(x, Pos, Strand, Score, Oligo, BG = FALSE,
        use_full_BG = FALSE, fullBG = FALSE, withDimnames = TRUE
    ) {
        x@ps_hits_pos <- Pos
        x@ps_hits_strand <- Strand
        x@ps_hits_score <- Score
        if (!use_full_BG) {
            x@ps_hits_score <- .ps_norm_score(x)
        }

        if (BG) {
            ps_bg_size(x) <- length(x@ps_hits_pos)
            ps_bg_avg(x) <- mean(x@ps_hits_score, na.rm = TRUE)
            ps_bg_std_dev(x) <- sd(x@ps_hits_score, na.rm = TRUE)
            if (ps_bg_std_dev(x) == 0) {
                ps_bg_std_dev(x) <- 0.00001
            }

            if (fullBG) {
                x@ps_hits_pos_bg <- Pos
                x@ps_hits_strand_bg <- Strand
                x@ps_hits_score_bg <- x@ps_hits_score
                names(x@ps_hits_score_bg) <- x@ps_bg_seq_names
                x@ps_hits_oligo_bg <- Oligo
            }
            x@ps_hits_pos <- integer()
            x@ps_hits_strand <- character()
            x@ps_hits_score <- numeric()
        } else {
            if (!is.na(x@ps_bg_avg) && !is.na(x@ps_bg_std_dev)) {
                # One-sample upper-tail z-test against the background mean.
                # The upper tail is evaluated directly with
                # pnorm(lower.tail = FALSE): computing it as 1 - pnorm(z)
                # saturates at 1.11e-16 and returns exactly 0 for z > ~8.3,
                # which silently discards the most enriched motifs.
                scores <- x@ps_hits_score[!is.na(x@ps_hits_score)]
                if (length(scores) <= 2L) {
                    stop("not enough foreground observations for the z-test",
                        call. = FALSE)
                }
                std_err <- x@ps_bg_std_dev / sqrt(length(scores))
                zscore <- (mean(scores) - x@ps_bg_avg) / std_err

                x@ps_zscore <- c(z = zscore)
                x@ps_pvalue <- pnorm(zscore, lower.tail = FALSE)
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
    }
)


setMethod(".ps_norm_score", "PSMatrix", function(x) {
    ps_score <- 1 + ((maxScore(Matrix(x)) - ps_hits_score(x)) /
        (minScore(Matrix(x)) - maxScore(Matrix(x))))

    return(ps_score)
})


setMethod(".ps_bg_from_table", "PSMatrix", function(x, short.matrix) {
    if (any(row.names(short.matrix) == ID(x))) {
        x@ps_bg_size <- as.integer(short.matrix[ID(x), "BG_SIZE"])
        x@ps_bg_avg <- as.numeric(short.matrix[ID(x), "BG_MEAN"])
        x@ps_bg_std_dev <- as.numeric(short.matrix[ID(x), "BG_STDEV"])
    } else {
        warning(sprintf("No background values found for %s %s", ID(x), name(x)))
    }

    return(x)
})


#' @importMethodsFrom TFBSTools Matrix
setMethod(".ps_norm_matrix", "PSMatrix", function(x) {
    mx <- Matrix(x)

    mx <- sweep(mx, 2, colSums(mx), FUN = "/")

    mx <- mx + x@.PS_PSEUDOCOUNT

    mx <- sweep(mx, 2, colSums(mx), FUN = "/")

    mx <- log(mx)

    Matrix(x) <- mx

    return(x)
})

# .ps_set_seq_names, .ps_scan_use_full_bg, and .ps_scan_standard are internal
# helpers used by ps_scan.
.ps_set_seq_names <- function(x, seqs, BG, use_full_BG) {
    if (BG || use_full_BG) {
        if (use_full_BG) {
            x@ps_seq_names <- seqs
        } else {
            x@ps_bg_seq_names <- names(seqs)
        }
    } else {
        x@ps_seq_names <- names(seqs)
    }
    x
}

.ps_scan_use_full_bg <- function(x, seqs, BG, use_full_BG) {
    # `seqs` are retained background sequence names resolved by pscan_fullBG(),
    # which are already exactly the names carried here. Matching them verbatim
    # avoids collapsing distinct transcripts that differ only after a dot.
    indices <- match(seqs, names(x@ps_hits_score_bg))
    res <- list(
        score = x@ps_hits_score_bg[indices],
        strand = x@ps_hits_strand_bg[indices],
        pos = x@ps_hits_pos_bg[indices],
        oligo = x@ps_hits_oligo_bg[indices]
    )
    x <- .ps_add_hits(x,
        Score = res$score, Strand = res$strand, Pos = res$pos,
        Oligo = res$oligo, BG = BG, use_full_BG = use_full_BG
    )
    x@ps_bg_seq_names <- character()
    x
}

# .ps_encode_seqs converts a character vector of EQUAL-LENGTH sequences into an
# (Nseq x L) integer matrix: A = 1, C = 2, G = 3, T = 4 (upper- or lower-case),
# every other character (e.g. N or an IUPAC code) -> NA. Row i corresponds to
# sequence i; encoded values match the A,C,G,T rows of the motif matrix.
# Returns NULL when the sequences are empty or not all the same length, which
# signals the caller to fall back to the per-sequence kernel .ps_scan_s.
# Doing the encoding once (in the multi-matrix callers) avoids re-encoding every
# sequence for every motif.
.ps_encode_bases <- function(seqs) {
    b <- as.integer(charToRaw(paste0(seqs, collapse = "")))
    code <- rep(NA_integer_, length(b))
    code[b == 65L | b == 97L] <- 1L
    code[b == 67L | b == 99L] <- 2L
    code[b == 71L | b == 103L] <- 3L
    code[b == 84L | b == 116L] <- 4L
    code
}

.ps_encode_seqs <- function(seqs) {
    n <- length(seqs)
    if (n == 0L) {
        return(NULL)
    }
    L <- nchar(seqs)
    if (L[1L] == 0L || any(L != L[1L])) {
        return(NULL)
    }
    code <- .ps_encode_bases(seqs)
    encoded <- matrix(code, nrow = n, ncol = L[1L], byrow = TRUE)
    has_ambiguous <- anyNA(code)
    attr(encoded, "has_ambiguous") <- has_ambiguous
    if (!has_ambiguous && length(code) <= 2e6) {
        # PWMscoreStartingAt scores every valid offset in one compiled call per
        # strand. Keep one subject with the shared encoding so that it is
        # built once for a multi-motif scan rather than once per motif. Large
        # background scans retain the bounded-memory R path below.
        attr(encoded, "subject") <- Biostrings::DNAString(
            paste0(seqs, collapse = "")
        )
    }
    encoded
}

# Score unambiguous, equal-length sequences through Biostrings' vectorised PWM
# primitive. Valid starts are selected within each concatenated sequence, so no
# window can cross a sequence boundary. Matrix layout and tie handling reproduce
# .ps_scan_s: first window on each strand, then forward strand on equal scores.
.ps_scan_biostrings <- function(seqs, subject, M, M_rc, W, n, L) {
    nw <- L - W + 1L
    subject_width <- n * L
    starting_at <- seq_len(subject_width - W + 1L)
    rownames(M) <- Biostrings::DNA_BASES
    rownames(M_rc) <- Biostrings::DNA_BASES

    score_fwd <- Biostrings::PWMscoreStartingAt(
        M, subject,
        starting.at = starting_at
    )
    score_rev <- Biostrings::PWMscoreStartingAt(
        M_rc, subject,
        starting.at = starting_at
    )

    # Pad the final positions so each sequence occupies one L-sized column.
    # Only the first nw rows are valid starts within a sequence.
    length(score_fwd) <- subject_width
    length(score_rev) <- subject_width
    dim(score_fwd) <- c(L, n)
    dim(score_rev) <- c(L, n)
    score_fwd <- score_fwd[seq_len(nw), , drop = FALSE]
    score_rev <- score_rev[seq_len(nw), , drop = FALSE]

    fwd_pos <- max.col(t(score_fwd), ties.method = "first")
    rev_pos <- max.col(t(score_rev), ties.method = "first")
    columns <- seq_len(n)
    fwd_value <- score_fwd[cbind(fwd_pos, columns)]
    rev_value <- score_rev[cbind(rev_pos, columns)]
    pick_fwd <- fwd_value >= rev_value

    best_pos <- rev_pos
    best_pos[pick_fwd] <- fwd_pos[pick_fwd]
    best_value <- rev_value
    best_value[pick_fwd] <- fwd_value[pick_fwd]
    best_strand <- rep("-", n)
    best_strand[pick_fwd] <- "+"

    list(
        score = best_value,
        strand = best_strand,
        pos = best_pos,
        oligo = unname(substring(seqs, best_pos, best_pos + W - 1L))
    )
}

# .ps_scan_batched scores every (equal-length) sequence against one motif in a
# vectorised way: for each of the W motif positions it adds a whole
# (sequence block x windows) slab of scores, then selects the best window per
# sequence. It reproduces .ps_scan_s exactly -- same column-sum order, same
# which.max() first-hit tie-breaking, same NA handling -- but without a separate
# R call per sequence. Sequences are processed in blocks to bound peak memory.
.ps_pick_block_hits <- function(score_fwd, score_rev, all_na) {
    n <- nrow(score_fwd)
    rows <- seq_len(n)
    fwd_pos <- max.col(score_fwd, ties.method = "first")
    rev_pos <- max.col(score_rev, ties.method = "first")
    fwd_value <- score_fwd[cbind(rows, fwd_pos)]
    rev_value <- score_rev[cbind(rows, rev_pos)]
    pick_fwd <- fwd_value >= rev_value

    pos <- rev_pos
    pos[pick_fwd] <- fwd_pos[pick_fwd]
    score <- rev_value
    score[pick_fwd] <- fwd_value[pick_fwd]
    strand <- rep("-", n)
    strand[pick_fwd] <- "+"
    if (any(all_na)) {
        score[all_na] <- NA_real_
        strand[all_na] <- "+"
        pos[all_na] <- 1L
    }
    list(score = score, strand = strand, pos = pos)
}

.ps_score_encoded_block <- function(S, M, M_rc, W, nw, has_ambiguous) {
    n <- nrow(S)
    score_fwd <- matrix(0, nrow = n, ncol = nw)
    score_rev <- matrix(0, nrow = n, ncol = nw)
    for (j in seq_len(W)) {
        columns <- seq.int(j, nw + j - 1L)
        bases <- S[, columns, drop = FALSE]
        score_fwd <- score_fwd + M[, j][bases]
        score_rev <- score_rev + M_rc[, j][bases]
    }
    if (!has_ambiguous) {
        return(.ps_pick_block_hits(score_fwd, score_rev, rep(FALSE, n)))
    }
    missing <- is.na(score_fwd)
    all_na <- rowSums(!missing) == 0L
    score_fwd[missing] <- -Inf
    score_rev[missing] <- -Inf
    .ps_pick_block_hits(score_fwd, score_rev, all_na)
}

.ps_scan_encoded_blocks <- function(seqs, S, M, M_rc, W, has_ambiguous) {
    n <- nrow(S)
    nw <- ncol(S) - W + 1L
    score <- numeric(n)
    pos <- integer(n)
    strand <- character(n)
    oligo <- character(n)
    block <- max(1L, as.integer(256e3 %/% nw))

    for (start in seq.int(1L, n, by = block)) {
        idx <- seq.int(start, min(start + block - 1L, n))
        hit <- .ps_score_encoded_block(
            S[idx, , drop = FALSE], M, M_rc, W, nw, has_ambiguous
        )
        score[idx] <- hit$score
        pos[idx] <- hit$pos
        strand[idx] <- hit$strand
        oligo[idx] <- substring(seqs[idx], hit$pos, hit$pos + W - 1L)
    }
    list(score = score, strand = strand, pos = pos, oligo = oligo)
}

.ps_scan_batched <- function(seqs, S, M, M_rc, W) {
    n <- nrow(S)
    L <- ncol(S)
    if (L < W) {
        return(list(
            score = rep(NA_real_, n), strand = rep("+", n),
            pos = rep(1L, n), oligo = rep("", n)
        ))
    }
    has_ambiguous <- attr(S, "has_ambiguous", exact = TRUE)
    subject <- attr(S, "subject", exact = TRUE)
    if (identical(has_ambiguous, FALSE) && !is.null(subject)) {
        return(.ps_scan_biostrings(seqs, subject, M, M_rc, W, n, L))
    }
    if (is.null(has_ambiguous)) {
        has_ambiguous <- anyNA(S)
    }
    .ps_scan_encoded_blocks(seqs, S, M, M_rc, W, has_ambiguous)
}

.ps_scan_standard <- function(
    x, seqs, BG, use_full_BG, fullBG, encoded = NULL
) {
    rc_x <- reverseComplement(x)
    nrows <- length(.PS_ALPHABET(x))
    W <- ncol(Matrix(x))
    # Forward and reverse-complement score matrices: built once per motif.
    M <- matrix(as.numeric(Matrix(x)), nrow = nrows, ncol = W)
    M_rc <- matrix(as.numeric(Matrix(rc_x)), nrow = nrows, ncol = W)

    if (is.null(encoded)) {
        encoded <- .ps_encode_seqs(seqs)
    }

    if (!is.null(encoded)) {
        # Fast path for equal-length sequences: vectorised batched scan.
        res <- .ps_scan_batched(seqs, encoded, M, M_rc, W)
    } else {
        # Fallback: sequences of differing length -> per-sequence kernel.
        Margs <- list(M = M, M_rc = M_rc, W = W)
        r <- mapply(.ps_scan_s, list(x), seqs, MoreArgs = Margs)
        res <- list(
            score = as.numeric(r["score", ]),
            strand = as.character(r["strand", ]),
            pos = as.integer(r["pos", ]),
            oligo = as.character(r["oligo", ])
        )
    }

    x <- .ps_add_hits(x,
        Score = res$score,
        Strand = res$strand,
        Pos = as.integer(res$pos),
        Oligo = res$oligo,
        BG = BG,
        use_full_BG = use_full_BG,
        fullBG = fullBG
    )
    if (!fullBG) {
        x@ps_bg_seq_names <- character()
    }
    x
}

#' Perform a Scan of DNA Sequences Using a `PSMatrix` Object
#'
#' This method performs a scan of DNA sequences using a `PSMatrix` object
#' to identify potential hits based on the matrix score. The scan is
#' performed on both strands (the direct and reverse complement strands) of the
#' sequences to ensure that no potential binding sites are missed.
#'
#' @param x A `PSMatrix` object representing the motif and related statistics
#'     to be used for scanning DNA sequences.
#' @param seqs A `DNAStringSet` object representing the sequences to be
#'    scanned for motif occurrences.
#' @param BG A logical value indicating whether to calculate background
#'    statistics.
#'    Default is set to `FALSE`.
#' @param use_full_BG A logical value (default is `FALSE`).
#'    If `TRUE`, the method assumes that the PSMatrix represents a special
#'    "full-background" set.
#'    In this case, the `seqs` parameter should be a vector of sequence names
#'    (instead of actual DNA sequences). The function handles this by matching
#'    sequence names against precomputed background data in the `PSMatrix`.
#' @param fullBG Logical, default is `FALSE`. If `TRUE` it computes a series
#'    of computation for the generation of a "full-background", a special case
#'    of background that retains all the background hits score, position,
#'    strand, and oligonucleotide sequence for each regulatory sequence scanned
#'    with a PWM.
#' @param encoded Internal. An optional pre-encoded integer matrix of the
#'    sequences (one column per sequence) produced by the calling functions to
#'    avoid re-encoding the same sequences for every motif. Users should leave
#'    this as the default `NULL`.
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
#' method matches the sequence names to those in the `PSMatrix`'s background
#' hit data and retrieves the corresponding binding information (score, strand,
#' position, oligo).
#'
#' When `use_full_BG` is set to `FALSE`, the function scans both the forward and
#' reverse complement strands of the sequences to ensure all potential binding
#' sites are detected. Optionally, the background statistics
#' (background average and standard deviation) can be computed and used during
#' scanning when `BG` is set to `TRUE`.
#'
#' Each motif position is scored against every window of the sequence; the
#' highest-scoring window on either strand is reported. Bases are matched
#' case-insensitively, so soft-masked (lower-case) sequence is scored as its
#' upper-case equivalent. Any window containing a character outside A/C/G/T
#' (for instance `N` or an IUPAC ambiguity code) is skipped. If a sequence is
#' shorter than the motif, or contains no scorable window at all (e.g. a fully
#' `N`-masked sequence), its score is reported as `NA` so that it does not
#' affect the foreground and background averages.
#'
#' @examples
#' matrix_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' matrix <- readRDS(matrix_path)
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)[1:5]
#'
#' scanned_result <- ps_scan(matrix, prom_seq)
#' scanned_result
#'
#' @export
#
#' @usage
#' \S4method{ps_scan}{PSMatrix}(x, seqs, BG = FALSE,
#'   use_full_BG = FALSE, fullBG = FALSE, encoded = NULL)
setMethod(
    "ps_scan",
    "PSMatrix",
    function(
        x, seqs, BG = FALSE, use_full_BG = FALSE, fullBG = FALSE,
        encoded = NULL
    ) {
        if (!is(seqs, "DNAStringSet") && !use_full_BG) {
            stop("seqs is not an object of DNAStringSet class")
        }

        x <- .ps_set_seq_names(x, seqs, BG, use_full_BG)
        seqs <- as.character(seqs)

        if (use_full_BG == TRUE) {
            x <- .ps_scan_use_full_bg(x, seqs, BG, use_full_BG)
        } else {
            x <- .ps_scan_standard(
                x, seqs, BG, use_full_BG, fullBG, encoded = encoded
            )
        }

        return(x)
    }
)


#' @importMethodsFrom Biostrings maxScore minScore reverseComplement
.ps_score_single_windows <- function(code, M, M_rc, W) {
    nw <- length(code) - W + 1L
    score_fwd <- numeric(nw)
    score_rev <- numeric(nw)
    for (j in seq_len(W)) {
        bases <- code[seq.int(j, nw + j - 1L)]
        score_fwd <- score_fwd + M[bases, j]
        score_rev <- score_rev + M_rc[bases, j]
    }
    list(forward = score_fwd, reverse = score_rev)
}

.ps_pick_single_hit <- function(scores, scores_rc, Seq, W) {
    fwd_pos <- which.max(scores)
    rev_pos <- which.max(scores_rc)
    if (length(fwd_pos) == 0L && length(rev_pos) == 0L) {
        return(list(
            score = NA_real_, strand = "+", pos = 1L,
            oligo = substring(Seq, 1L, W)
        ))
    }
    if (scores[fwd_pos] >= scores_rc[rev_pos]) {
        return(list(
            score = scores[fwd_pos], strand = "+", pos = fwd_pos,
            oligo = substring(Seq, fwd_pos, fwd_pos + W - 1L)
        ))
    }
    list(
        score = scores_rc[rev_pos], strand = "-", pos = rev_pos,
        oligo = substring(Seq, rev_pos, rev_pos + W - 1L)
    )
}

setMethod(".ps_scan_s", "PSMatrix", function(x, Seq, M, M_rc, W) {
    if (nchar(Seq) < W) {
        return(list(score = NA_real_, strand = "+", pos = 1L, oligo = ""))
    }
    scores <- .ps_score_single_windows(.ps_encode_bases(Seq), M, M_rc, W)
    .ps_pick_single_hit(scores$forward, scores$reverse, Seq, W)
})

# setMethod(".ps_assign_score", "PSMatrix", function(x, S){
# sum(Matrix(x)[matrix(data = c(.PS_ALPHABET(x)[S], 1:length(x)), ncol = 2, nrow
# = length(x))])
# })

# .ps_assign_score scored a single window one base at a time. It is no longer
# used: .ps_scan_s now scores all windows at once via matrix indexing. Kept
# (commented out) for reference.
# .ps_assign_score <- function(S, x, AB, ncolx) {
#   sum(x[ncolx + AB[S]]) # Assign score to oligo
# }

#' Validate a PSMatrix object
#'
#' This function checks if the input `PSMatrix` object properties meet the ones
#' desired.
#'
#' @param object A `PSMatrix` object.
#'
#' @details
#' The function `validPSMatrix` ensures that the `PSMatrix` object has
#' correctly formatted background and foreground averages, standard deviations,
#' and hit-related values.
#' It performs the following steps:
#' \itemize{
#'   \item `ps_bg_avg`,`ps_fg_avg`, and `ps_bg_std_dev` must be of length 1.
#'   \item The values of \code{ps_bg_avg} and \code{ps_bg_std_dev}
#'   must be between 0 and 1 (excluding 0 for \code{ps_bg_std_dev}).
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
validPSMatrix <- function(object) {
    if (length(object@ps_bg_avg) != 1) {
        return("Background average must be of length 1")
    }
    if (length(object@ps_fg_avg) != 1) {
        return("Foreground average must be of length 1")
    }
    if (length(object@ps_bg_std_dev) != 1) {
        return("Background stdev must be of length 1")
    }
    if ((object@ps_bg_avg < 0 || object@ps_bg_avg > 1) &&
        !is.na(object@ps_bg_avg)) {
        return(paste(
            "Invalid value for Background average: ", object@ps_bg_avg
        ))
    }
    if ((object@ps_bg_std_dev < 0 || object@ps_bg_std_dev > 1) &&
        !is.na(object@ps_bg_std_dev)) {
        return(paste(
            "Invalid value for Background stddev: ", object@ps_bg_std_dev
        ))
    }
    # if(object@ps_bg_size < 1000 && !is.na(object@ps_bg_size))
    # return(paste("Invalid value for Background size: ", object@ps_bg_size, "
    # Background must be of at least 1000 sequences"))
    if (length(object@ps_hits_pos) != length(object@ps_hits_strand) ||
        length(object@ps_hits_pos) != length(object@ps_hits_score)) {
        return(paste(
            "Invalid PSMatrix object: different values for hits, strands ",
            "and scores vectors"
        ))
    }

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
#'    \code{\link{ps_fg_size}}, \code{\link{ps_zscore}},
#'    \code{\link{ps_pvalue}}
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
#' This method specifically sets the background average for an object of
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
#' This method allows users to modify this value.
#'
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_avg<-`(pfm1, value = 0.73473)
#'
#' @export
setReplaceMethod("ps_bg_avg", "PSMatrix", function(x, value) {
    x@ps_bg_avg <- value
    validObject(x)
    x
})

#' Set Background Standard Deviation
#'
#' This method updates the background standard deviation stored in a `PSMatrix`
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
setReplaceMethod("ps_bg_std_dev", "PSMatrix", function(x, value) {
    x@ps_bg_std_dev <- value
    validObject(x)
    x
})

#' Set Background Size in a PSMatrix
#'
#' This method specifically sets the background size for a `PSMatrix` object.
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
setReplaceMethod("ps_bg_size", "PSMatrix", function(x, value) {
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
setAs("PFMatrix", "PSMatrix", function(from) {
    # .ps_norm_matrix(new(
    #   "PSMatrix", from, ps_bg_avg = as.numeric(NA),
    #   ps_fg_avg = as.numeric(NA), ps_bg_std_dev = as.numeric(NA),
    #   ps_bg_size = as.integer(NA), .PS_PSEUDOCOUNT = 0.01
    # ))

    PSMatrix(from)
})

# PSMatrix <- function(pfm, ps_bg_avg = as.numeric(NA), ps_fg_avg =
# as.numeric(NA), ps_bg_std_dev = as.numeric(NA),
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
setAs("PFMatrixList", "PSMatrixList", function(from) {
    to <- lapply(from, as, "PSMatrix")

    # A PFMatrixList has no legend, but a PSMatrixList coerced through here
    # does, and dropping it would lose a full background at the first step of
    # ps_build_bg().
    do.call(PSMatrixList, c(to, list(transcriptIDLegend = .ps_legend_of(from))))
})
