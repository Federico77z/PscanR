#' Transcript ID Legend (Generic Function)
#'
#' `transcriptIDLegend` is a **generic function** that retrieves the
#' `transcriptIDLegend` slot from a `PSMatrixList` object. This slot contains
#' a character vector mapping transcript IDs used in the dataset.
#'
#' @param x A `PSMatrixList` object.
#' @param ... Additional arguments (not used in this method but included
#'    for extensibility).
#'
#' @return A character vector containing transcript ID mappings.
#'
#' @export
#'
#' @examples
#' full_pfms_path <- system.file("extdata", "full_pfms.rds", package = "PscanR")
#' full_pfms <- readRDS(full_pfms_path)
#' head(transcriptIDLegend(full_pfms))
#'
setGeneric(
    "transcriptIDLegend",
    function(x, ...) standardGeneric("transcriptIDLegend")
)

#' Retrieve the Motif-Enrichment Z-Score (Generic Function)
#'
#' `ps_zscore` retrieves the Z-statistic stored in a `PSMatrix` after scanning.
#'
#' The Z-score represents how significantly a motif is enriched in
#' the input sequences compared to the background model.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result): a **numeric value**
#'   representing the stored Z-score.
#' }
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_zscore(pfm1)
#'
#' @export
setGeneric("ps_zscore", function(x, ...) standardGeneric("ps_zscore"))

#' Retrieve the Motif-Enrichment P-Value (Generic Function)
#'
#' `ps_pvalue` retrieves the upper-tail p-value stored in a `PSMatrix` after
#' scanning.
#'
#' The p-value quantifies the statistical significance of motif
#' enrichment in the input sequences compared to the background model.
#' A lower p-value is stronger evidence against the background null model.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result): a **numeric value**
#'   representing the stored p-value.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_pvalue(pfm1)
#'
#' @export
setGeneric("ps_pvalue", function(x, ...) standardGeneric("ps_pvalue"))

#' Retrieve the Background Average Score for Motif Enrichment Analysis
#' (Generic Function)
#'
#' `ps_bg_avg` retrieves the average background score stored in a `PSMatrix`.
#'
#' The background average score represents the mean motif score computed
#' over the background sequences. This value is essential for normalizing
#' motif enrichment results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result): a **numeric value**
#'   representing the stored background average score.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_avg(pfm1)
#'
#' @export
setGeneric("ps_bg_avg", function(x, ...) standardGeneric("ps_bg_avg"))

#' Retrieve the Foreground Average Score for Motif Enrichment Analysis
#' (Generic Function)
#'
#' `ps_fg_avg` retrieves the average foreground score stored in a `PSMatrix`.
#'
#' The foreground average score represents the mean motif score computed
#' over the input (foreground) sequences. This value is crucial for
#' comparing motif occurrences between the foreground and background models.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result): a **numeric value**
#'   representing the stored foreground average score.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_fg_avg(pfm1)
#'
#' @export
setGeneric("ps_fg_avg", function(x, ...) standardGeneric("ps_fg_avg"))

#' Retrieve the Background Standard Deviation for Motif Enrichment
#' Analysis (Generic Function)
#'
#' `ps_bg_std_dev` retrieves the background standard deviation stored in a
#' `PSMatrix`.
#'
#' The background standard deviation value quantifies the variability in the
#' background motif binding scores.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result): a **numeric value**
#'   representing the stored background standard deviation value.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_std_dev(pfm1)
#'
#' @export
setGeneric("ps_bg_std_dev", function(x, ...) standardGeneric("ps_bg_std_dev"))

#' Retrieve the Background Size Value for Motif Enrichment Analysis
#' (Generic Function)
#'
#' `ps_bg_size` is a **generic function** that retrieves the background
#' size value stored in a `PSMatrix`.
#'
#' The background size is the number of sequences in the background reference
#' set used to calculate the stored statistics.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result), an **integer value**.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_bg_size(pfm1)
#'
#' @export
setGeneric("ps_bg_size", function(x, ...) standardGeneric("ps_bg_size"))

#' Retrieve the Foreground Size Value for Motif Enrichment Analysis
#' (Generic Function)
#'
#' `ps_fg_size` is a **generic function** that retrieves the foreground
#' size value stored in a `PSMatrix`.
#'
#' The foreground size is the number of unique sequences retained after input
#' cleaning, including sequences whose best-hit score is `NA`.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result), an **integer value**.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_fg_size(pfm1)
#'
#' @export
setGeneric("ps_fg_size", function(x, ...) standardGeneric("ps_fg_size"))

#' Retrieve the Hits Size Value for Motif Enrichment Analysis
#' (Generic Function)
#'
#' `ps_hits_size` is a **generic function** that retrieves the number of hits
#' for motif enrichment analysis. Methods should be
#' implemented for specific object classes that store scan results.
#'
#' The hits size is the number of retained per-sequence best windows. It does
#' not apply a score or significance threshold.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object (Pscan result), an **integer value**.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_size(pfm1)
#'
#' @export
setGeneric("ps_hits_size", function(x, ...) standardGeneric("ps_hits_size"))

#' Retrieve Motif Hit Scores (Generic Function)
#'
#' `ps_hits_score` is a **generic function** that retrieves the motif hit
#' scores for each promoter sequence in an object. These scores represent
#' the binding affinity or enrichment level of promoter sequences when
#' scanned with a Position Weight Matrix (PWM).
#'
#' Methods should be implemented for specific object classes that store
#' motif scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **named numeric vector**, where
#'   names correspond to promoter sequence identifiers and values represent
#'   their respective motif hit scores.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_score(pfm1)
#'
#' @export
setGeneric("ps_hits_score", function(x, ...) standardGeneric("ps_hits_score"))

#' Retrieve Motif Hit Scores from a Background Dataset (Generic Function)
#'
#' ps_hits_score_bg is a generic function designed to extract motif hit scores
#' from a background dataset. These scores quantify the binding affinity or
#' enrichment level of promoter sequences when scanned using a Position Weight
#' Matrix (PWM). This function is specifically useful for analyzing background
#' matrices, which serve as reference datasets in pscan analyses.
#'
#' Methods should be implemented for specific object classes that store
#' motif scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' The `ps_hits_score_bg` function is particularly relevant for background
#' datasets, which include stored motif scan results for every retained
#' background sequence. These scores provide a reference for comparison in
#' motif enrichment analyses, helping to assess the significance of observed
#' motif occurrences.
#'
#' In a PSMatrixList object, the `ps_hits_score_bg` slot is populated only for
#' full-background matrices. This allows `pscan_fullBG()` to retrieve selected
#' hits without rescanning DNA.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **named numeric vector**, where
#'   names correspond to promoter sequence identifiers and values represent
#'   their respective motif hit scores.
#' }
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'   "full_pfm1.rds",
#'   package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_score_bg(full_pfm1)
#'
#' @export
setGeneric(
    "ps_hits_score_bg",
    function(x, ...) standardGeneric("ps_hits_score_bg")
)

#' Compute Motif Hit Z-Scores (Generic Function)
#'
#' `ps_hits_z` is a **generic function** that computes Z-scores for motif hit
#' scores in an object. The Z-score indicates how unusual a motif score is
#' compared with scores in the supplied background reference set.
#' Higher Z-scores suggest stronger motif enrichment, which may indicate
#' regulatory significance.
#'
#' Methods should be implemented for specific object classes that store
#' motif scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **named numeric vector**, where
#'   names correspond to promoter sequence identifiers and values represent
#'   their respective Z-scores.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_z(pfm1)
#'
#' @export
setGeneric("ps_hits_z", function(x, ...) standardGeneric("ps_hits_z"))

#' Retrieve Motif Hit Strand Information (Generic Function)
#'
#' `ps_hits_strand` is a **generic function** that retrieves the strand
#' information (`+` or `-`) of motif hits within an object. This indicates
#' whether a motif was detected on the forward (`+`) or reverse (`-`) strand
#' of a promoter sequence.
#'
#' The Pscan algorithm scans both strands of promoter sequences to ensure that
#' no potential transcription factor binding sites are missed.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector**, where
#'   names correspond to promoter sequence identifiers, and values represent
#'   the strand (`+` or `-`) on which the motif was detected.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_strand(pfm1)
#'
#' @export
setGeneric(
    "ps_hits_strand",
    function(x, ...) standardGeneric("ps_hits_strand")
)

#' Retrieve Motif Hit Strand Information for the Background Dataset
#' (Generic Function)
#'
#' `ps_hits_strand_bg` is a **generic function** that retrieves the strand
#' information (`+` or `-`) of motif hits within an object. This indicates
#' whether a motif was detected on the forward (`+`) or reverse (`-`) strand
#' of a promoter sequence of the background dataset.
#'
#' The Pscan algorithm scans both strands of promoter sequences to ensure that
#' no potential transcription factor binding sites are missed.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' The `ps_hits_strand_bg` function is particularly relevant for background
#' datasets, which include stored motif scan results for every retained
#' background sequence.
#'
#' In a PSMatrixList object, the `ps_hits_strand_bg` slot is populated only
#' for full-background matrices. This allows `pscan_fullBG()` to retrieve
#' selected hits without rescanning DNA.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector**, where
#'   names correspond to promoter sequence identifiers, and values represent
#'   the strand (`+` or `-`) on which the motif was detected.
#' }
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'   "full_pfm1.rds",
#'   package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_strand_bg(full_pfm1)
#'
#' @export
setGeneric(
    "ps_hits_strand_bg",
    function(x, ...) standardGeneric("ps_hits_strand_bg")
)

#' Retrieve Motif Hit Positions (Generic Function)
#'
#' `ps_hits_pos` is a **generic function** that retrieves the positions of
#' motif hits in a given object. These positions indicate where motifs are
#' located within promoter sequences.
#' For the `PSMatrix` method, the optional `pos_shift` is the coordinate of the
#' first sequence base; stored one-based positions are converted accordingly.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: an **integer vector**, where
#'   names correspond to promoter sequence identifiers, and values represent
#'   motif hit positions (with an optional shift).
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_pos(pfm1)
#'
#' @export
setGeneric("ps_hits_pos", function(x, ...) standardGeneric("ps_hits_pos"))

#' Retrieve Motif Hit Positions on Background Dataset (Generic Function)
#'
#' `ps_hits_pos_bg` is a **generic function** that retrieves the positions of
#' motif hits in a given object. These positions indicate where motifs are
#' located within promoter sequences of a background dataset.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' The `ps_hits_pos_bg` function is particularly relevant for background
#' datasets, which include stored motif scan results for every retained
#' background sequence.
#'
#' In a PSMatrixList object, the `ps_hits_pos_bg` slot is populated only for
#' full-background matrices. This allows `pscan_fullBG()` to retrieve selected
#' hits without rescanning DNA.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: an **integer vector**, where
#'   names correspond to promoter sequence identifiers, and values represent
#'   motif hit positions along the background dataset.
#' }
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'   "full_pfm1.rds",
#'   package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_pos_bg(full_pfm1)
#'
#' @export
setGeneric("ps_hits_pos_bg", function(x, ...) standardGeneric("ps_hits_pos_bg"))

#' Retrieve Matched Oligonucleotide Sequences (Generic Function)
#'
#' `ps_hits_oligo` is a **generic function** that extracts the oligonucleotide
#' sequences corresponding to motif hits in a given object.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector**, where
#'   names correspond to sequence identifiers, and values represent the
#'   oligonucleotide sequences (subset of the input promoter sequences)
#'   matching the motif.
#' }
#'
#' @section Strand:
#' The sequence returned is always the **forward strand** of the promoter, as
#' it appears in the input, whatever strand the motif matched on. Reading it
#' against the matrix consensus is therefore misleading for a reverse-strand
#' hit: pair it with \code{\link{ps_hits_strand}} and reverse-complement the
#' oligo where that reports \code{"-"}.
#'
#' This is the behaviour of the original Pscan implementation and is kept for
#' consistency with it. It is also the more useful convention for sequence
#' context, since oligos from different motifs are then all in the same frame
#' of reference and can be compared to each other and to the promoter.
#'
#' @seealso \code{\link{ps_hits_strand}}, \code{\link{ps_hits_pos}},
#'   \code{\link{ps_hits_table}}
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_oligo(pfm1)
#'
#' # Oligos are forward-strand, so reverse-strand hits must be flipped before
#' # they are compared with the motif consensus.
#' oligos <- ps_hits_oligo(pfm1)
#' reverse <- ps_hits_strand(pfm1) == "-"
#' oligos[reverse] <- as.character(Biostrings::reverseComplement(
#'     Biostrings::DNAStringSet(oligos[reverse])
#' ))
#' head(oligos)
#'
#' @export
setGeneric("ps_hits_oligo", function(x, ...) standardGeneric("ps_hits_oligo"))

#' Retrieve Matched Oligonucleotide Sequences from a Background Dataset
#' (Generic Function)
#'
#' `ps_hits_oligo_bg` is a **generic function** that extracts the
#' oligonucleotide sequences corresponding to motif hits among a background
#' dataset in a given object.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' The `ps_hits_oligo_bg` function is particularly relevant for background
#' datasets, which include stored motif scan results for every retained
#' background sequence.
#'
#' In a PSMatrixList object, the `ps_hits_oligo_bg` slot is populated only
#' for full-background matrices. This allows `pscan_fullBG()` to retrieve
#' selected hits without rescanning DNA.
#'
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector**, where
#'   names correspond to sequence identifiers, and values represent the
#'   oligonucleotide sequences (subset of the promoter sequences) matching
#'   the motif.
#' }
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'   "full_pfm1.rds",
#'   package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_hits_oligo_bg(full_pfm1)
#'
#' @export
setGeneric(
    "ps_hits_oligo_bg",
    function(x, ...) standardGeneric("ps_hits_oligo_bg")
)

#' Retrieve a Summary Table of Motif Hits (Generic Function)
#'
#' `ps_hits_table` is a **generic function** that extracts motif hit
#' information from an object and organizes it into a structured table.
#' The table includes motif hit scores, strand orientation, positions,
#' and the matched oligonucleotide sequences.
#'
#' Methods should be implemented for specific object classes that store motif
#' scanning results.
#'
#' @param x An object containing motif scan results.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' - If `x` is a `PSMatrix` object: a **data frame** with the following
#'   columns:
#'   \itemize{
#'     \item `SCORE`: The motif hit score.
#'     \item `POS`: The position of the motif hit.
#'     \item `STRAND`: The strand orientation (`+` or `-`).
#'     \item `OLIGO`: The matched oligonucleotide sequence.
#'   }
#'   Rows correspond to sequence names and are **sorted by decreasing score**.
#'
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_hits_table(pfm1)
#'
#' @export
setGeneric("ps_hits_table", function(x, ...) standardGeneric("ps_hits_table"))

#' Retrieve Sequence Names (Generic Function)
#'
#' `ps_seq_names` is a **generic function** that extracts the sequence names
#' or identifiers from an object containing promoter sequence data.
#'
#' Methods should be implemented for specific object classes storing sequence
#' data.
#'
#' @param x An object containing sequence information.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector** of sequence
#'   names corresponding to the analyzed promoter regions.
#' }
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' ps_seq_names(pfm1)
#'
#' @export
setGeneric("ps_seq_names", function(x, ...) standardGeneric("ps_seq_names"))

#' Retrieve Sequence Names for the Background Dataset (Generic Function)
#'
#' `ps_bg_seq_names` is a **generic function** that extracts the sequence names
#' or identifiers stored for a full-background reference set.
#'
#' Methods should be implemented for specific object classes storing sequence
#' data.
#'
#' @param x An object containing sequence information.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' \itemize{
#'   \item If `x` is a `PSMatrix` object: a **character vector** of sequence
#'   names corresponding to the analyzed promoter regions.
#' }
#'
#' @details
#' The `ps_bg_seq_names` function is particularly relevant for background
#' datasets, which include stored motif scan results for every retained
#' background sequence.
#'
#' In a PSMatrixList object, the `ps_bg_seq_names` slot is populated only for
#' full-background matrices. This allows `pscan_fullBG()` to retrieve selected
#' hits without rescanning DNA.
#'
#' @examples
#' full_pfm1_path <- system.file("extdata",
#'   "full_pfm1.rds",
#'   package = "PscanR"
#' )
#' full_pfm1 <- readRDS(full_pfm1_path)
#' ps_bg_seq_names(full_pfm1)
#'
#' @export
setGeneric(
    "ps_bg_seq_names",
    function(x, ...) standardGeneric("ps_bg_seq_names")
)

setGeneric(
    ".PS_PSEUDOCOUNT",
    function(x, ...) standardGeneric(".PS_PSEUDOCOUNT")
)


setGeneric(
    ".PS_ALPHABET",
    function(x, ...) standardGeneric(".PS_ALPHABET")
)


setGeneric(
    ".ps_norm_matrix",
    function(x, ...) standardGeneric(".ps_norm_matrix")
)

setGeneric(
    ".ps_bg_seq_names",
    function(x, out) standardGeneric(".ps_bg_seq_names")
)

setGeneric(".ps_seq_names", function(x, out) standardGeneric(".ps_seq_names"))

#' Perform Motif Scanning (Generic Function)
#'
#' `ps_scan` is a **generic function** that scans DNA sequences for motif
#' occurrences using an object containing motif scoring data. The scan is
#' typically performed on both the forward and reverse complement strands to
#' detect all potential transcription factor binding sites.
#'
#' Methods should be implemented for specific object classes that support motif
#' scanning.
#'
#' @param x An object containing motif information, typically a `PSMatrix`
#'    object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return
#' - If `x` is a `PSMatrix` object: a **PSMatrix object** with updated motif
#'   hit information, including positions, strands, scores, and matched
#'   sequences.
#'
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
setGeneric("ps_scan", function(x, ...) standardGeneric("ps_scan"))


setGeneric(
    ".ps_bg_from_table",
    function(x, ...) standardGeneric(".ps_bg_from_table")
)


setGeneric(".ps_scan_s", function(x, ...) standardGeneric(".ps_scan_s"))


setGeneric(
    ".ps_norm_score",
    function(x, ...) standardGeneric(".ps_norm_score")
)

# setGeneric(".ps_assign_score",
#           function(x, ...) standardGeneric(".ps_assign_score"))


setGeneric(
    ".ps_add_hits",
    function(x, ...) standardGeneric(".ps_add_hits")
)

setGeneric(
    "all_sequences_ID",
    function(x, ...) standardGeneric("all_sequences_ID")
)

#' Generic Setter Method for Background Mean Value of an object
#'
#' @param x An object for which the background average is to be set.
#' @param ... Additional arguments.
#' @param value The value for the background average.
#'
#' @return An object of the same class of the input with the modified
#'    background average value.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_avg<-`(pfm1, value = 0.73473)
#'
#' @export
setGeneric(
    "ps_bg_avg<-",
    function(x, ..., value) standardGeneric("ps_bg_avg<-")
)

#' Set Background Standard Deviation (Generic Function)
#'
#' `ps_bg_std_dev<-` is a **generic function** that sets the background
#' standard deviation value for a given object.
#'
#' Methods should be implemented for specific object classes that support
#' background standard deviation assignment.
#'
#' @param x An object to which the background standard deviation will be
#'    assigned, typically a `PSMatrix` object.
#' @param ... Additional arguments (not used in the generic function but may be
#'    used in specific methods).
#' @param value A numeric value representing the background standard deviation.
#'
#' @return The updated object with the assigned background standard deviation.
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_std_dev<-`(pfm1, value = 0.08568925)
#'
#' @export
setGeneric(
    "ps_bg_std_dev<-",
    function(x, ..., value) standardGeneric("ps_bg_std_dev<-")
)

#' Set Background Size (Generic Function)
#'
#' `ps_bg_size<-` is a **generic function** that assigns the background
#' size for a given object. The background size typically represents
#' the number of background sequences used for motif enrichment analysis.
#'
#' Methods should be implemented for specific object classes that support
#' background size assignment.
#'
#' @param x An object to which the background size will be assigned,
#'    typically a `PSMatrix` object.
#' @param ... Additional arguments (not used in the generic function but may be
#'    used in specific methods).
#' @param value A numeric value representing the background size (e.g., number
#'    of background sequences).
#'
#' @return The updated object with the assigned background size.
#'
#' @examples
#' pfm1_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' pfm1 <- readRDS(pfm1_path)
#' `ps_bg_size<-`(pfm1, value = 25629L)
#'
#' @export
setGeneric(
    "ps_bg_size<-",
    function(x, ..., value) standardGeneric("ps_bg_size<-")
)
