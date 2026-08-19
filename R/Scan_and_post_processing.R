#' Executes the Pscan algorithm on a set of regulatory sequences.
#'
#' This function computes alignment scores between regulatory sequences
#' (a set of gene promoters, `x`) and position weight matrices,
#' which quantify potential binding affinities of transcription factors in that
#' region.
#' These matrices can be sourced from public databases such as JASPAR.
#' The scanning is performed by the Pscan algorithm.
#'
#' @param x A `DNAStringSet` object containing the set of regulatory sequences
#'    from co-regulated or co-expressed genes (i.e. a set of gene promoters).
#'    See the Biostrings package for details.
#'
#' @param pfms A `PSMatrixList` object containing PWMs and background
#'    statistics.
#'    For background statistics we refer to the standard deviation and average
#'    of hits scores when the background (set of promoters of all the
#'    transcript in the organism of study) is scanned with the position weight
#'    matrices.
#'    This is used to assert the statistical enrichment of motif occurrences
#'    in co-expressed or co-regulated genes.
#'    See \code{\link{ps_build_bg}}, \code{\link{ps_retrieve_bg_from_file}},
#'    \code{\link{ps_build_bg_from_table}},
#'    \code{\link{generate_psmatrixlist_from_background}} for how to create
#'    `PSMatrixList` objects that contain background statistics.
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
#'    processing many tasks. log = TRUE enables logging to debug each step of
#'    the parallel tasks.
#'
#' @details
#' The `pscan` function performs sequence scanning using the `ps_scan` method
#' for individual PWMs, accounting for both the forward and reverse complement
#' strands, ensuring no potential binding sites are missed.
#' The method extracts all k-mer substrings of the input sequences with the
#' same length of the motif and evaluates a score for both the forward and
#' reverse strand based on how well the k-mer matches the Transcription Factor
#' Binding Motif provided by the PWM. For each input sequence scanned with
#' individual PWM, only the k-mer with the highest score is selected.
#' The method populates the following PSMatrix slots for each input
#' sequence:
#' \itemize{
#'    \item ps_hits_score: the highest matching value found for each sequence.
#'    \item ps_hits_pos: position of the best TFBS match.
#'    \item ps_hits_strand: the DNA strand of the TFBS ('+' for forward and '-'
#'    for reverse).
#'    \item ps_hits_oligo: the sequence of the binding site.}
#'
#' This function uses example datasets located in the `extdata/` directory for
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the
#' examples.
#'
#' @return
#' Enriched PSMatrixList object where each matrix includes
#' \code{ps_hits_score}, \code{ps_hits_pos}, \code{ps_hits_strand},
#' and \code{ps_hits_oligo} populated for each input sequence.
#'
#' @export
#'
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions:
#' # -200 +50 bp in respect to the TSS.
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a small background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' motif_ids <- c(
#'   "MA0506.1", "MA0632.2", "MA0615.1",
#'   "MA0076.2", "MA0645.1", "MA0838.1"
#' )
#' background <- background[motif_ids]
#'
#' # Execute the PScan algorithm
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_results_table(results)
#'
pscan <- function(x, pfms, BPPARAM = bpparam(), BPOPTIONS = bpoptions()) {
    .ps_checks(x, pfms, type = 4)

    x <- BiocGenerics::unique(x)
    x <- .clean_sequence(x)

    # Encode the sequences once and reuse the encoding for every motif.
    encoded <- .ps_encode_seqs(as.character(x))

    pfms <- BiocParallel::bplapply(
    pfms,
    FUN = ps_scan,
    x,
    BG = FALSE,
    encoded = encoded,
    BPPARAM = BPPARAM,
    BPOPTIONS = BPOPTIONS
    )

    BiocGenerics::do.call(PSMatrixList, pfms)
}

#' Extract Pre-Computed Metrics from a PSMatrixList Object
#'
#' This function is a faster alternative to \code{\link{pscan}}, that works
#' only if a full background PSMatrixList object has already been computed.
#' Instead of re-scanning input sequences, it retrieves pre-computed alignment
#' scores, TFBS positions, strands, and binding oligos for a selected subset of
#' sequences from the complete background.
#'
#' This is particularly useful when analyzing different subsets of sequences
#' on the same background, reducing computation time.
#'
#' @param ID A character vector containing sequence identifiers (transcript
#'    IDs) corresponding to the sequences of interest. These should match
#'    those used when the background was originally computed
#'    (e.g., RefSeq transcript IDs).
#' @param full_pfms A `PSMatrixList` object generated using
#'    \code{\link{ps_build_bg}} with the `fullBG = TRUE` flag. This object
#'    must contain a complete set of pre-computed scores and metadata for all
#'    promoter sequences of the reference organism.
#' @param scheme Transcript identifier scheme. `"auto"` (the default) detects
#'    it from the background's transcript identifiers. See
#'    \code{\link{ps_select_promoters}} for the available schemes.
#' @param quiet Logical. Suppress the scheme-detection message.
#'
#' @details
#' Internally, the function matches the provided `ID` values against the
#' `transcriptIDLegend` slot of the full background PSMatrixList.
#' The transcriptIDLegend slot was introduced to handle cases where duplicate
#' promoter sequences are removed using unique() during background generation.
#' Since in this case only one representative sequence name is retained,
#' the mapping vector ensures that if a user provides a transcript ID
#' corresponding to a removed duplicate, it can still be correctly associated
#' with the retained sequence used in the background.
#'
#' Matching is scheme-aware. Under the RefSeq scheme a trailing version
#' (`".1"`, `".2"`) is a release of the same transcript and is normalised away,
#' so `NM_000546.6` matches a background built from `NM_000546.4`. Under the
#' TAIR scheme the same punctuation denotes a splice variant, so
#' `AT1G01110.2` is a different transcript from `AT1G01110.1` and the suffix is
#' preserved. Requesting an identifier that is ambiguous within the background
#' raises a warning naming the identifier rather than silently returning the
#' first match.
#'
#' Identifiers absent from the background are dropped with a warning. That
#' normally means the promoter was excluded during background construction
#' (high 'N' content or a length mismatch), but it also covers identifiers that
#' were never part of this background.
#'
#' @return A `PSMatrixList` object in which the alignment scores and related
#'    metrics have been retrieved from `full_pfms` for each sequence in `ID`,
#'    based on the background data.
#'
#' @examples
#' # load a full background PSMatrixList
#' full_bg_path <- system.file("extdata", "full_pfms.rds", package = "PscanR")
#' full_pfms <- readRDS(full_bg_path)
#'
#' IDs <- c(
#'   "NM_031921.6", "NM_001005484.2", "NR_047525.1", "NR_029639.1",
#'   "NR_036051.1", "NR_029834.1", "NM_001029885.2", "NR_148357.1",
#'   "NR_148960.1", "NM_001130413.4"
#' )
#'
#' results <- pscan_fullBG(IDs, full_pfms)
#' ps_results_table(results)
#'
#' @export
pscan_fullBG <- function(ID, full_pfms, scheme = "auto", quiet = FALSE) {
    if (!is.character(ID)) {
    stop("ID must be a character vector containing transcript identifiers")
    }
    if (length(full_pfms@transcriptIDLegend) == 0) {
    stop("The background PSMatrixList must be a full background")
    }

    all_seq_ID <- full_pfms@transcriptIDLegend

    # Reduce both sides to the transcript key the scheme defines. For RefSeq
    # that removes a version; for TAIR it removes nothing, because there the
    # suffix distinguishes splice variants with different TSSs. The legend
    # *values* are left untouched: they already name sequences in the
    # background exactly.
    scheme <- .ps_resolve_scheme(
        scheme = scheme, promoter_ids = names(all_seq_ID),
        annotation_ids = ID, quiet = quiet
    )
    legend_key <- .ps_transcript_base(names(all_seq_ID), scheme)
    query_key <- .ps_transcript_base(ID, scheme)

    # A key that names more than one background transcript cannot be resolved.
    # Taking the first silently is how the wrong promoter used to be returned.
    duplicated_keys <- unique(legend_key[duplicated(legend_key)])
    ambiguous <- intersect(query_key, duplicated_keys)
    if (length(ambiguous) > 0) {
    warning(
        "Ambiguous transcript identifier(s) under scheme \"", scheme,
        "\": ", paste(utils::head(ambiguous, 10), collapse = ", "),
        if (length(ambiguous) > 10) ", ..." else "",
        ". Each matches more than one background transcript; the first was ",
        "used. Supply fully qualified identifiers to disambiguate.",
        call. = FALSE
    )
    }

    # Use of all_seq_ID (mapping vector) to extract sequences name retained
    # in full BG that have the same sequence to those inserted by the user.
    x <- stats::setNames(unname(all_seq_ID)[match(query_key, legend_key)], ID)

    # NA removal
    rem_names <- names(x[is.na(x)])

    if (length(rem_names) > 0) {
    warning(paste(
        "Found", length(rem_names), "identifier(s) absent from the background",
        "(excluded during background construction for high N content or a",
        "length mismatch, or never part of it). Removing:",
        paste(rem_names, collapse = ", ")
    ))
    }

    x <- x[!is.na(x)]

    .check_seq_duplicated(x)

    x <- unique(x)

    # See ps_scan for details

    pfms <- lapply(
    full_pfms,
    FUN = ps_scan,
    x,
    BG = FALSE,
    use_full_BG = TRUE
    )

    BiocGenerics::do.call(PSMatrixList, pfms)
}

# .ps_check_filtered_inputs and .ps_filter_promoters are internal helpers used
# by PscanFiltered.
.ps_check_filtered_inputs <- function(prom_seq, Jmatrix, background) {
    if (!is(prom_seq, "DNAStringSet")) {
    stop("Invalid input: 'prom_seq' must be a DNAStringSet")
    }
    if (!is(Jmatrix, "PSMatrix")) {
    stop("Invalid input: 'Jmatrix' must be a PSMatrix")
    }
    if (!is(background, "PSMatrixList")) {
    stop("Invalid input: 'background' must be a PSMatrixList")
    }
}

.ps_filter_promoters <- function(prom_seq, Jmatrix, n) {
    threshold <- ps_bg_avg(Jmatrix) + n * ps_bg_std_dev(Jmatrix)
    rc_matrix <- reverseComplement(Jmatrix)
    nrows <- length(.PS_ALPHABET(Jmatrix))
    W <- ncol(Matrix(Jmatrix))
    # Build the strand score matrices once and reuse them for every sequence.
    M <- matrix(as.numeric(Matrix(Jmatrix)), nrow = nrows, ncol = W)
    M_rc <- matrix(as.numeric(Matrix(rc_matrix)), nrow = nrows, ncol = W)
    seqs <- as.character(prom_seq)
    encoded <- .ps_encode_seqs(seqs)
    if (!is.null(encoded)) {
    Jmatrix@ps_hits_score <- .ps_scan_batched(seqs, encoded, M, M_rc, W)$score
    } else {
    res <- mapply(.ps_scan_s, list(Jmatrix), seqs,
        MoreArgs = list(M = M, M_rc = M_rc, W = W)
    )
    Jmatrix@ps_hits_score <- as.numeric(res["score", ])
    }
    Jmatrix@ps_hits_score <- .ps_norm_score(Jmatrix)
    Score <- Jmatrix@ps_hits_score

    filtered_prom_seq <- character()
    if (n > 0) {
    filtered_prom_seq <- prom_seq[Score >= threshold]
    }
    if (n < 0) {
    filtered_prom_seq <- prom_seq[Score <= threshold]
    }
    if (length(filtered_prom_seq) == 0) {
    warning("No sequence satisfy the filter criterium")
    return(NULL)
    }
    filtered_prom_seq
}

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
#' @param background A PSMatrixList object. It contains PWMs and background
#'    statistics.
#'    For background statistics we refer to the standard deviation and average
#'    of hits scores when the background (set of promoters of all the
#'    transcript in the organism of study) is scanned with the position weight
#'    matrices.
#'    This is used to assert the statistical enrichment of motif occurrences
#'    in co-expressed or co-regulated genes.
#'    See \code{\link{ps_build_bg}}, \code{\link{ps_retrieve_bg_from_file}},
#'    \code{\link{ps_build_bg_from_table}},
#'    \code{\link{generate_psmatrixlist_from_background}}
#'    for how to create `PSMatrixList` objects that contain background
#'    statistics.
#' @param BPPARAM The BPPARAM used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`.
#'    If BPPARAM is not explicitly set, the default value (bpparam()) will be
#'    used, which automatically chooses a sensible parallelization method based
#'    on the user's system.
#'    You can specify BPPARAM = BiocParallel::SnowParam(8) on all operating
#'    systems, or BPPARAM = BiocParallel::MulticoreParam(8) on Unix-like
#'    systems to use, for example, 8 cores.
#' @param BPOPTIONS The BPOPTIONS used by bplapply. See BiocParallel package.
#'    This argument is passed to `BiocParallel::bplapply`.
#'    The default is `bpoptions()`.
#'    Some useful tasks: bpoptions(progressbar = TRUE, log = TRUE).
#'    progressbar = TRUE enables a progress bar that can be useful when
#'    processing many tasks. log = TRUE enables logging to debug each step of
#'    the parallel tasks.
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
#'
#' @examples
#' # Run on a small subset using bundled data (no downloads).
#' prom_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(prom_path)
#' set.seed(1)
#' n_prom <- min(100, length(prom_seq))
#' prom_seq <- prom_seq[sample(seq_along(prom_seq), n_prom)]
#'
#' bg_path <- system.file("extdata",
#'   "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#' bg <- ps_retrieve_bg_from_file(bg_path, J2020)
#'
#' JM <- bg[[1]]
#' res <- PscanFiltered(prom_seq,
#'   JM,
#'   background = bg,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' @export
PscanFiltered <- function(prom_seq, Jmatrix, n = 1, background,
                            BPPARAM = bpparam(), BPOPTIONS = bpoptions()) {
    .ps_check_filtered_inputs(prom_seq, Jmatrix, background)

    filtered_prom_seq <- .ps_filter_promoters(prom_seq, Jmatrix, n)
    if (is.null(filtered_prom_seq)) {
    return(NULL)
    }

    pfms <- BiocParallel::bplapply(
    background,
    FUN = ps_scan,
    filtered_prom_seq,
    BG = FALSE,
    BPPARAM = BPPARAM,
    BPOPTIONS = BPOPTIONS
    )

    BiocGenerics::do.call(PSMatrixList, pfms)
}

#' Create a Summary Table of PscanR Results
#'
#' This function generates a table summarizing the statistical analysis of
#' motif enrichment, stored in a `PSMatrixList` object. It retrieves
#' pre-computed background and foreground statistics (the alignment scores
#' between Pscan input regulatory sequences and position weight matrices)
#' for each PWM.
#' Then it compares motif occurrences in the foreground (e.g.,
#' coexpressed/coregulated promoter sequences) with the background (e.g., all
#' promoters in an organism) using Z-scores, p-values, and FDR correction.
#' Returns a table ordered by decreasing `ZSCORE` and increasing `P.VALUE`.
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated
#'    metadata (foreground and background statistics). Typically the output
#'    of `pscan()` or `pscan_fullBG()` functions.
#' @param FDR A numeric value indicating the maximum false discovery
#'    rate (FDR) allowed for filtering the result table.
#'    Only rows with FDR <= `FDR` will be retained.
#'    Default is 1, so all the results are displayed.
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
#'   (FG_AVG-BG_AVG)/BG_STDEV. It represents a statistical measure of motif
#'   enrichment, quantifying how much the average foreground score deviates from
#'   the background one.
#'   \item "P.VALUE": The p-value for each PWM.
#'   \item "FDR": The adjusted p-value, representing the False Discovery Rate,
#'   is calculated using the Benjamini-Hochberg correction.
#' }
#'
#' @details
#'
#' This summary table is useful for the identification of transcription factors
#' that could be common regulators of a set of co-expressed or co-regulated
#' genes.
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
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a small background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' motif_ids <- c(
#'   "MA0506.1", "MA0632.2", "MA0615.1",
#'   "MA0076.2", "MA0645.1", "MA0838.1"
#' )
#' background <- background[motif_ids]
#'
#' # Execute the PScan algorithm and view the result table
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_results_table(results, FDR = 0.1)
#'
#' @export
#' @importFrom stats p.adjust
ps_results_table <- function(pfms, FDR = 1) {
    .ps_checks2(pfms)

    if (!is.numeric(FDR) || length(FDR) != 1L || is.na(FDR) ||
        FDR < 0 || FDR > 1) {
    stop("FDR must be a single numeric value between 0 and 1")
    }

    bg_v <- vapply(pfms, ps_bg_avg, numeric(length = 1L))
    std_v <- vapply(pfms, ps_bg_std_dev, numeric(length = 1L))
    fg_v <- vapply(pfms, ps_fg_avg, numeric(length = 1L))
    zs_v <- vapply(pfms, ps_zscore, numeric(length = 1L))
    pv_v <- vapply(pfms, ps_pvalue, numeric(length = 1L))
    fdr_v <- p.adjust(pv_v, method = "BH")

    tbl <- data.frame(
    "NAME" = name(pfms), "BG_AVG" = bg_v, "BG_STDEV" = std_v,
    "FG_AVG" = fg_v, "ZSCORE" = zs_v,
    "P.VALUE" = pv_v, "FDR" = fdr_v, row.names = ID(pfms)
    )

    tbl <- tbl[tbl$FDR <= FDR, , drop = FALSE]

    tbl[with(tbl, order(P.VALUE, ZSCORE, decreasing = c(FALSE, TRUE))),
        , drop = FALSE]
}

#' Generate a Z-score Matrix of Motif Enrichment per Sequence
#'
#' This function creates a matrix of Z-scores for motif occurrences across
#' multiple PSMatrix contained in a `PSMatrixList` object. Each row corresponds
#' to a sequence, and each column represents a motif.
#'
#' @param pfms A `PSMatrixList` object containing multiple Position Weight
#'    Matrices and associated metadata (foreground and background statistics).
#'    This object is the output of `pscan()` function.
#'
#' @details
#' Z-score represents the statistical significance of the alignment scores
#' for regulatory sequences relative to background expectation. A high Z-score
#' suggests strong motif enrichment in the foreground compared to
#' the background.
#'
#' This function uses example datasets located in the `extdata/` directory for
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the
#' examples.
#'
#' @return
#' A matrix in which each column corresponds to a motif in the `PSMatrixList`
#' object, and each row to the sequence identifiers. The matrix contains
#' z-score values for each sequence-motif pair, indicating how much the
#' computed scores deviate from the expected background distribution.
#' Higher Z-score values represent stronger motif enrichment.
#'
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions:
#' # -200 +50 bp in respect to the TSS.
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a small background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' motif_ids <- c(
#'   "MA0506.1", "MA0632.2", "MA0615.1",
#'   "MA0076.2", "MA0645.1", "MA0838.1"
#' )
#' background <- background[motif_ids]
#'
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_z_table(results)
#'
#' @export
ps_z_table <- function(pfms) {
    .ps_checks2(pfms)

    tbl <- lapply(pfms, ps_hits_z)

    # as.matrix(as.data.frame(tbl, col.names = name(pfms)))

    as.matrix(as.data.frame(tbl, col.names = ID(pfms)))
}

#' Motif Z-Score Correlation Heatmap from Pscan Results
#'
#' This function generates a heatmap that visualizes the correlation of
#' Z-scores for transcription factors (TFs) based on a specified false
#' discovery rate threshold.
#' It allows customization of the heatmap's appearance.
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated
#'    metadata (foreground and background statistics). Typically the output
#'    of `pscan()` or `pscan_fullBG()` functions.
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to include in the analysis. The default is `0.01`.
#' @param ... Additional user-defined arguments to customize the heatmap
#'    settings, such as color palettes or clustering objects.
#'
#' @details
#' The heatmap represents the correlation of motif Z-scores, helping to
#' identify clusters of motifs that show similar enrichment patterns across
#' sequences.
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts the result table and the z-score table from the pscan
#'      algorithm result.
#'   \item Filters the result table based on the given FDR threshold.
#'   \item Uses the filtered result table to create the z-table for the
#'      selected TFs.
#'   \item Generates the heatmap using the `pheatmap` function,
#'      with customizable settings.}
#'
#'  Default settings (which can be changed) are:
#'  \itemize{
#'    \item `cluster_rows` and `cluster_cols`, set to `TRUE` by default.
#'    \item `color`, which uses a blue-white-red palette.
#'    \item The main `Pscan Score Correlation Heatmap`.
#'    \item TF names are shown as column labels.
#'    }
#'
#' This function uses example datasets located in the `extdata/` directory for
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the
#' examples.
#'
#' @return Invisibly returns the Z-score submatrix used for plotting.
#'   The heatmap is drawn as a side effect.
#'
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions:
#' # -200 +50 bp in respect to the TSS.
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a small background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' motif_ids <- c(
#'   "MA0506.1", "MA0632.2", "MA0615.1",
#'   "MA0076.2", "MA0645.1", "MA0838.1"
#' )
#' background <- background[motif_ids]
#'
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_zscore_heatmap(results, FDR = 0.05)
#'
#' @export
#' @import pheatmap
#' @importFrom utils modifyList
#' @importFrom grDevices colorRampPalette
ps_zscore_heatmap <- function(pfms, FDR = 0.01, ...) {
    res_table <- ps_results_table(pfms)
    z_table <- ps_z_table(pfms)
    topn <- which(res_table$FDR <= FDR)

    defaults <- list(
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = "Pscan Score Correlation Heatmap", scale = "row",
    show_rownames = FALSE,
    labels_col = res_table$NAME[topn],
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    fontsize = 10,
    fontsize_row = 6,
    fontsize_col = 6
    )

    user_args <- list(...)

    final_args <- modifyList(defaults, user_args)

    tf_to_plot <- rownames(res_table)[topn]

    z_table_reduced <- z_table[, tf_to_plot]

    res <- do.call(pheatmap::pheatmap, c(list(z_table_reduced), final_args))

    invisible(z_table_reduced)
}

#' Heatmap of Motif Hit Positions from Pscan Results
#'
#' This function creates a heatmap visualizing the positional distribution
#' of motif hits based on a specified false discovery rate threshold.
#' The heatmap helps visualize where significant motif hits occur within the
#' analyzed sequences.
#'
#' @param pfms A `PSMatrixList` object containing multiple PWMs and associated
#'    metadata (foreground and background statistics). Typically the output
#'    of `pscan()` or `pscan_fullBG()` function.
#' @param FDR Numeric. False Discovery Rate (FDR) threshold to select the TFs
#'    to be included in the analysis. The default is set to `0.01`.
#' @param shift Integer. A value to shift the reported positions of motif hits
#'    in respect to the TSS.
#'    Default is set to `0`.
#' @param ... Additional user-defined arguments that can be passed to
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
#' Default settings (which can be changed) are:
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
#' @return Invisibly returns the positional hits matrix used for plotting.
#'   The heatmap is drawn as a side effect.
#'
#' @examples
#' # Load the promoter sequences for hg38 (Homo sapiens), promoter regions:
#' # -200 +50 bp in respect to the TSS.
#' file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a small background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' motif_ids <- c(
#'   "MA0506.1", "MA0632.2", "MA0615.1",
#'   "MA0076.2", "MA0645.1", "MA0838.1"
#' )
#' background <- background[motif_ids]
#'
#' # Execute the Pscan algorithm and view the result table
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_hitpos_map(results, shift = -200)
#'
#' @export
#' @import pheatmap
#' @importFrom utils modifyList
#' @importFrom grDevices colorRampPalette
ps_hitpos_map <- function(pfms, FDR = 0.01, shift = 0, ...) {
    res_table <- ps_results_table(pfms)

    topn <- which(res_table$FDR <= FDR)

    defaults <- list(
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "yellow", "red"))(100),
    main = "Pscan Hits Position Heatmap",
    fontsize = 10, show_rownames = FALSE, scale = "none",
    clustering_distance_rows = "manhattan",
    clustering_distance_cols = "manhattan",
    clustering_method = "average",
    fontsize_row = 6,
    fontsize_col = 6
    )

    user_args <- list(...)

    final_args <- modifyList(defaults, user_args)

    pos_mat <- matrix(
    data = NA, nrow = ps_fg_size(pfms[[1]]),
    ncol = length(topn)
    )

    # `topn` indexes rows of res_table; the matrix columns are 1..length(topn).
    # Conflating the two only works while `topn` happens to be 1..k, which is
    # true today because BH-adjusted FDR is monotone in the p-value ordering
    # ps_results_table() sorts by, but it is not a property to rely on.
    for (i in seq_along(topn)) {
    pos_mat[, i] <- ps_hits_pos(pfms[[row.names(res_table)[topn[i]]]],
        pos_shift = shift
    )
    }

    colnames(pos_mat) <- res_table$NAME[topn]
    rownames(pos_mat) <- ps_seq_names(pfms[[1]])

    res <- do.call(pheatmap::pheatmap, c(list(pos_mat), final_args))

    invisible(pos_mat)
}

# Shared styling for the ggplot2 plotters, so that figures produced by
# different functions in this file read as one family.
#
# Categorical fill colours, chosen so that adjacent pairs stay distinguishable
# under the common forms of colour-vision deficiency. Beyond this many groups
# ps_motif_barplot() uses the ggplot2 default instead, since no fixed order
# stays separable.
.PS_GROUP_PALETTE <- c(
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#4a3aa7", "#e34948"
)

# Single-series colour for the density figures, taken from the same palette.
.PS_DENSITY_COLOUR <- .PS_GROUP_PALETTE[[1L]]

# The three score bands the two background thresholds cut a scan into, weakest
# first. Named in the vocabulary the `st` argument uses elsewhere, so a reader
# can carry the legend straight across to ps_density_plot().
.PS_SCORE_BANDS <- c(
    "below background mean",
    "above mean (loose)",
    "above mean + SD (strict)"
)

#' Density Plot of Motif Hits Along Promoter Regions
#'
#' This function creates a density plot representing the distribution of hits
#' along the promoter sequences based on their position and score for a
#' specific PWM.
#' It helps visualize where motif occurrences are concentrated.
#'
#' @param pfm A `PSMatrix` object. It is a selected matrix from the
#'    `PSMatrixList`, result of `pscan` function.
#' @param shift Integer value specifying the positional shift applied to the
#'    hit positions to obtain the position in respect to the TSS.
#'    Default is `0`.
#' @param st Score threshold used to filter hits. Can be a numeric value to set
#'    the threshold directly, or a character:
#'    \itemize{
#'      \item `all`: the threshold is set to `0` (All the hits are evaluated).
#'      \item `loose`: uses the background average score as threshold.
#'      \item `strict`: uses the background average score together with the
#'      background standard deviation as threshold.}
#'    Default is `loose` (background average).
#' @param window Optional length-2 numeric giving the promoter window, on the
#'    same scale as `shift` -- for promoters extracted 1000 bp upstream of the
#'    TSS and plotted with `shift = -1000`, that is `c(-1000, 0)`. Supplying it
#'    corrects the estimate at the edges by reflection; leaving it `NULL` only
#'    clips the curve to the observed range. The reflecting bounds are the
#'    positions a hit can be *reported* at, which stop `ncol(pfm) - 1` short of
#'    the end of the window, since a hit is recorded at the start of the motif.
#' @param bins Optional number of bins for a binned profile of the same
#'    positions, drawn behind the curve. Bars are on the density scale, so they
#'    share the curve's y axis directly. Carrying no bandwidth and no boundary
#'    assumption, they show what the smoothing is doing.
#'
#' @return A `ggplot` object showing the density of hit positions along the
#'    promoter region. Nothing is drawn until the object is printed, so it can
#'    be modified with further `ggplot2` layers first.
#'
#' @details
#' The function filters motif hits based on a specified threshold and generates
#' a density plot to show their distribution. The function includes a vertical
#' dashed line marking the mode (the most frequent position along the
#' promoters). The title reports how many promoters passed the threshold.
#'
#' This function uses example datasets located in the `extdata/` directory for
#' demonstration purposes only. These files are not part of the core data used
#' by the function. They can be accessed using `system.file()` as shown in the
#' examples.
#'
#' @examples
#' matrix_path <- system.file("extdata", "pfm1.rds", package = "PscanR")
#' matrix <- readRDS(matrix_path)
#' ps_density_plot(matrix, shift = -200)
#'
#' # The result is a ggplot object, so it can be restyled before drawing.
#' ps_density_plot(matrix, shift = -200) +
#'     ggplot2::labs(subtitle = "-200 to +50 around the TSS")
#'
#' @export
#' @importFrom stats density
ps_density_plot <- function(pfm, shift = 0, st = ps_bg_avg(pfm),
                            window = NULL, bins = NULL) {
    # st = score threshold. It can be passed as a numeric value
    # or as one of three characters "all", "loose", "strict".
    st <- .ps_resolve_threshold(st, pfm, "st")

    scores <- ps_hits_score(pfm)
    g_scores <- scores >= st
    sum_g <- sum(g_scores)

    positions <- ps_hits_pos(pfm, pos_shift = shift)[g_scores]
    # A hit is reported at the start of its window, so the last position a
    # motif of width w can start at is w - 1 short of the end of the promoter.
    # Reflecting at the end of the promoter itself would place the boundary
    # where no hit can be observed and flatten the estimate short of it.
    density_hits <- .ps_bounded_density(
    positions, .ps_hit_support(window, ncol(pfm))
    )

    .ps_density_ggplot(
    density_hits,
    title = paste(
        name(pfm), "binding site density across", sum_g, "promoter regions"
    ),
    xlab = "Position along promoters",
    values = positions, bins = bins
    )
}

.ps_resolve_threshold <- function(st, M, label) {
    if (!is.character(st)) {
    return(st)
    }
    switch(st,
    "all" = 0,
    "loose" = ps_bg_avg(M),
    "strict" = ps_bg_avg(M) + ps_bg_std_dev(M),
    {
        warning(sprintf("Invalid value for %s, reverting to loose", label))
        ps_bg_avg(M)
    }
    )
}

#' Bin values into a profile on the density scale
#'
#' Heights are counts divided by the total and by the bin width, which is the
#' same scale a kernel density is on. Both can then share one y axis, with no
#' secondary axis and no arbitrary rescaling, and the comparison between them is
#' meaningful rather than merely visual.
#'
#' The bars carry no bandwidth and no boundary assumption, so they also serve as
#' a check on the smoothed curve: where the two disagree, it is the curve that
#' is imposing something.
#'
#' @param values Numeric vector, or `NULL`.
#' @param bins Number of bins, or `NULL`.
#' @param limits Length-2 numeric spanning the bins.
#'
#' @return A data frame of `mid`, `density` and `width`, or `NULL`.
#'
#' @noRd
.ps_binned_profile <- function(values, bins, limits) {
    if (is.null(bins) || is.null(values)) {
    return(NULL)
    }
    if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || bins < 1) {
    stop("bins must be a single positive number", call. = FALSE)
    }
    breaks <- seq(limits[[1L]], limits[[2L]], length.out = as.integer(bins) + 1L)
    width <- diff(breaks)[[1L]]
    counts <- tabulate(
    cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE),
    nbins = as.integer(bins)
    )
    data.frame(
    mid = breaks[-length(breaks)] + width / 2,
    density = counts / (sum(counts) * width),
    width = width
    )
}

#' Positions a hit of a given width can be reported at within a window
#'
#' A hit is recorded at the start of the matching window, so a motif of width
#' `w` can never be reported in the final `w - 1` bases of the promoter.
#' Reflecting at the end of the promoter itself would put the boundary where no
#' hit can be observed.
#'
#' @param window Promoter window, or `NULL`.
#' @param width Motif width, `ncol()` of the matrix.
#'
#' @return Length-2 numeric, or `NULL` when `window` is `NULL`.
#'
#' @noRd
.ps_hit_support <- function(window, width) {
    if (is.null(window)) {
    return(NULL)
    }
    if (length(window) != 2L || !all(is.finite(window)) ||
        window[[1L]] == window[[2L]]) {
    stop(
        "'window' must be two finite, distinct values giving the promoter ",
        "window the positions were taken from",
        call. = FALSE
    )
    }
    window <- sort(as.numeric(window))
    support <- c(window[[1L]] + 1, window[[2L]] - width + 1)
    if (support[[1L]] >= support[[2L]]) {
    stop(
        "'window' spans ", diff(window), " bases, which is too short to ",
        "report a hit for a motif of width ", width,
        call. = FALSE
    )
    }
    support
}

#' Estimate a density on a bounded interval
#'
#' \code{stats::density()} assumes the quantity is unbounded. It evaluates over
#' \code{cut = 3} bandwidths beyond the data, so the curve runs outside the
#' interval the values can occupy, and it lets kernel mass escape across the
#' boundary, so the estimate is depressed at both ends. For hit positions the
#' interval is the scanned window, and both artefacts land exactly where
#' positional signal is usually read.
#'
#' With `window` supplied the estimate is corrected by reflection: the sample is
#' mirrored across both bounds before estimation, which returns the escaped mass
#' to the interval and makes the result integrate to one over it.
#'
#' Reflection is only meaningful at a real boundary. Without `window` the
#' support is unknown -- the sample minimum and maximum are not boundaries, they
#' are just the smallest and largest values seen -- so no correction is applied
#' and the grid is merely clipped to the observed range.
#'
#' Reflection constrains the estimate to have zero derivative at each bound. On
#' a distribution still changing at the edge this flattens the last bandwidth or
#' so, which is visible when the bandwidth is large relative to the structure.
#'
#' @param x Numeric vector.
#' @param window Optional length-2 numeric giving the interval `x` is confined
#'   to.
#'
#' @return A \code{stats::density} object.
#'
#' @noRd
#' @importFrom stats density
.ps_bounded_density <- function(x, window = NULL) {
    if (is.null(window)) {
    limits <- range(x)
    if (!all(is.finite(limits)) || limits[[1L]] == limits[[2L]]) {
        return(density(x))
    }
    return(density(x, from = limits[[1L]], to = limits[[2L]]))
    }

    window <- sort(as.numeric(window))
    if (length(window) != 2L || !all(is.finite(window)) ||
        window[[1L]] == window[[2L]]) {
    stop(
        "'window' must be two finite, distinct values giving the interval ",
        "the positions are confined to",
        call. = FALSE
    )
    }
    outside <- x < window[[1L]] | x > window[[2L]]
    if (any(outside)) {
    stop(
        sum(outside), " value(s) fall outside 'window' [",
        window[[1L]], ", ", window[[2L]], "]",
        call. = FALSE
    )
    }

    # Estimate on the sample mirrored across both bounds, then keep the middle
    # copy's worth of density. bw is taken from the original sample, since the
    # mirrored one is three times the size and would be over-smoothed.
    bw <- density(x)$bw
    mirrored <- c(x, 2 * window[[1L]] - x, 2 * window[[2L]] - x)
    out <- density(
    mirrored,
    bw = bw, from = window[[1L]], to = window[[2L]]
    )
    out$y <- out$y * 3
    out$n <- length(x)
    out
}

#' Draw a density object as a filled curve with its mode marked
#'
#' Both density plotters draw exactly this figure and differ only in their
#' labels. The curve is built from the evaluation grid of an existing
#' \code{stats::density} object rather than with \code{geom_density()}, so the
#' plotted curve is the one already computed -- same bandwidth, same grid -- and
#' the mode is read off the same object.
#'
#' @param d A \code{stats::density} object.
#' @param title,xlab Plot title and x axis label.
#' @param values The values `d` was estimated from, needed only for `bins`.
#' @param bins Number of bins for a binned profile drawn behind the curve, or
#'   `NULL` for none.
#'
#' @return A `ggplot` object.
#'
#' @noRd
#' @importFrom ggplot2 .data
.ps_density_ggplot <- function(d, title, xlab, values = NULL, bins = NULL) {
    curve <- data.frame(x = d$x, y = d$y)
    peak <- d$x[which.max(d$y)]
    bars <- .ps_binned_profile(values, bins, range(d$x))

    # Keep the mode label inside the panel: the base-graphics original always
    # wrote it to the right of the line, which clipped whenever the mode fell
    # near the end of the range.
    past_middle <- peak > mean(range(d$x))

    plot <- ggplot2::ggplot(curve, ggplot2::aes(x = .data$x, y = .data$y))
    plot <- if (is.null(bars)) {
    # Without bars the filled area gives the curve some weight.
    plot + ggplot2::geom_area(fill = .PS_DENSITY_COLOUR, alpha = 0.15)
    } else {
    # With bars it would only muddy them, so the curve is drawn on its own.
    plot + ggplot2::geom_col(
        data = bars,
        mapping = ggplot2::aes(x = .data$mid, y = .data$density),
        width = bars$width[[1L]], fill = "grey82", colour = "white",
        linewidth = 0.25
    )
    }

    plot +
    ggplot2::geom_line(colour = .PS_DENSITY_COLOUR, linewidth = 0.8) +
    ggplot2::geom_vline(
        xintercept = peak, colour = "grey50", linetype = "dashed",
        linewidth = 0.6
    ) +
    ggplot2::annotate(
        "text",
        x = peak, y = max(d$y), label = paste("Mode:", round(peak)),
        hjust = if (past_middle) 1.1 else -0.1, vjust = -0.6, size = 3.5,
        colour = "grey30"
    ) +
    # The label sits just above the peak, so the panel needs headroom that the
    # default expansion does not leave.
    ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(title = title, x = xlab, y = "Density") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 11)
    )
}

.ps_hit_distances <- function(M1, M2, st1, st2) {
    scores1 <- ps_hits_score(M1)
    g_scores1 <- scores1 >= st1
    scores2 <- ps_hits_score(M2)
    g_scores2 <- scores2 >= st2

    hits1 <- ps_hits_pos(M1, pos_shift = ncol(M1) / 2)[g_scores1]
    hits2 <- ps_hits_pos(M2, pos_shift = ncol(M2) / 2)[g_scores2]

    hits1 <- hits1[names(hits1) %in% names(hits2)]
    hits2 <- hits2[names(hits2) %in% names(hits1)]

    hits1 - hits2
}

#' Density Plot of Distances between Identified Motif Hits in two PSMatrix
#' Objects
#'
#' This function visualizes the density plot of distances between identified
#' hit sites in two `PSMatrix` objects. The distance between hits is calculated
#' for each sequence that is present in both matrices.
#' It allows filtering the identified sites based on a specified threshold
#' value.
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
#'      \item `strict`: uses the background average plus the background
#'      standard deviation as threshold.}
#'    Default is set to `loose`.
#' @param st2 Score threshold used to filter hits for the second `PSMatrix`.
#'    Can be a numeric value to set the threshold directly, or a character, as
#'    for st1.
#'    Default is set to `loose`.
#' @param window Optional length-2 numeric giving the interval the distances
#'    can occupy. Supplying it corrects the estimate at the edges by
#'    reflection; leaving it `NULL` only clips the curve to the observed range.
#' @param bins Optional number of bins for a binned profile of the same
#'    distances, drawn behind the curve on the density scale.
#'
#' @return A `ggplot` object showing the distribution of distances between
#'    identified motif hits in \code{M1} and \code{M2}. Nothing is drawn until
#'    the object is printed, so it can be modified with further `ggplot2`
#'    layers first.
#'
#' @seealso \code{\link{ps_bg_avg}}, \code{\link{ps_bg_std_dev}},
#' \code{\link{ps_hits_score}}, \code{\link{ps_hits_pos}}
#'
#' @details
#' The x-axis represents the distances between corresponding hits, computed as
#' the position of the hit in M1 minus the position of the hit in M2. Positions
#' run along the promoter in the direction of transcription, so a positive
#' value places M1 at the larger coordinate: downstream of M2, on the side of
#' the transcription start site. A negative value places M1 upstream of M2.
#' The y-axis represents the density of those distances, and a dashed line
#' marks the mode.
#'
#' Hit positions are shifted by half of the motif length (i.e., `ncol(M1)/2`
#' and `ncol(M2)/2`) before distances are computed, so distances are centered
#' on the motif midpoint.
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
#' prom_seq <- readRDS(file_path)[1:10]
#'
#' # Build a two-motif background from bundled files.
#' matrix_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' matrices <- readRDS(matrix_path)
#' bg_path <- system.file(
#'   "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'   package = "PscanR"
#' )
#' background <- ps_retrieve_bg_from_file(bg_path, matrices)
#' background <- background[c("MA0506.1", "MA0632.2")]
#'
#' results <- pscan(prom_seq, background,
#'   BPPARAM = BiocParallel::SerialParam()
#' )
#'
#' ps_density_distances_plot(results[[1]], results[[2]], "all", "loose")
#' @export
ps_density_distances_plot <- function(M1, M2, st1 = ps_bg_avg(M1),
                                        st2 = ps_bg_avg(M2),
                                        window = NULL, bins = NULL) {
    if (!is(M1, "PSMatrix") || !is(M2, "PSMatrix")) {
    stop("Both object must be of class PSMatrix")
    }

    st1 <- .ps_resolve_threshold(st1, M1, "st1")
    st2 <- .ps_resolve_threshold(st2, M2, "st2")
    distances <- .ps_hit_distances(M1, M2, st1, st2)
    density_distances <- .ps_bounded_density(distances, window)

    # The count is the point of the title as much as the names are: a distance
    # needs a retained hit for both motifs, so raising either threshold can
    # leave the plot resting on far fewer promoters than were scanned.
    .ps_density_ggplot(
    density_distances,
    title = paste(
        M1@name, "binding site distance relative to", M2@name,
        "across", length(distances), "promoter regions"
    ),
    xlab = "Distances between the identified sites",
    values = distances, bins = bins
    )
}

#' Structural Class of Each Motif
#'
#' Returns the structural class of the transcription factor behind each motif,
#' for example `"Basic leucine zipper factors (bZIP)"` or `"WRKY"`.
#'
#' @param pfms A `PSMatrixList` object, or any `PFMatrixList`.
#'
#' @return A character vector of class labels, named by JASPAR matrix
#'    identifier. Motifs with no class in JASPAR are reported as
#'    `"Unclassified"`.
#'
#' @details
#' The class is read from the `matrixClass` slot with
#' \code{TFBSTools::matrixClass()}. It is worth knowing where it lives, because
#' `tags(x)$class` is `NULL` for every JASPAR matrix in every release:
#' `TFBSTools` moves the class out of the tag list and into its own slot when
#' the matrix is built. Since `PSMatrix` extends `PFMatrix`, the slot survives
#' scanning, so the class is available from a scan result without consulting
#' JASPAR again.
#'
#' A small number of JASPAR matrices carry two classes; only the first is
#' reported, so the result is always one label per motif and is safe to use as
#' a grouping variable. Surrounding whitespace is stripped, because JASPAR
#' ships a handful of classes with a trailing space — in JASPAR 2020 the IRF3,
#' IRF4 and IRF5 matrices carry `"Tryptophan cluster factors "`, which would
#' otherwise group and colour separately from the other 38 members of that
#' class.
#'
#' Class is preferred over the `family` tag for grouping. Family coverage is
#' incomplete in some collections, and its vocabulary changes between JASPAR
#' releases, so it does not give a stable grouping across analyses.
#'
#' @seealso \code{\link{ps_motif_barplot}}
#'
#' @export
#'
#' @examples
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#' head(ps_motif_class(J2020))
ps_motif_class <- function(pfms) {
    if (!is(pfms, "PFMatrixList")) {
    stop("pfms must be a PSMatrixList or PFMatrixList object", call. = FALSE)
    }

    out <- vapply(pfms, function(x) {
    value <- matrixClass(x)
    if (length(value) == 0L || is.na(value[[1L]])) {
        return("Unclassified")
    }
    # JASPAR ships a few classes with stray surrounding whitespace, which
    # would otherwise split one class into two grouping levels.
    label <- trimws(value[[1L]])
    if (!nzchar(label)) {
        return("Unclassified")
    }
    label
    }, character(1L))

    stats::setNames(out, ID(pfms))
}

#' Bar Plot of the Most Enriched Motifs
#'
#' Plots the highest ranking motifs of a scan as horizontal bars, optionally
#' coloured by a grouping such as the structural class of the factor.
#'
#' @param pfms A `PSMatrixList` returned by \code{\link{pscan}}, or a
#'    `data.frame` produced by \code{\link{ps_results_table}}. A results table
#'    carries no motif metadata, so `group` must be supplied explicitly when
#'    one is used.
#' @param n Number of motifs to show. If the collection holds fewer than `n`
#'    motifs, all of them are shown.
#' @param statistic Column to rank and plot. One of `"ZSCORE"` (the default),
#'    `"P.VALUE"`, `"FDR"`, `"FG_AVG"` or `"BG_AVG"`.
#' @param group Optional grouping mapped to bar colour. `NULL` (the default)
#'    draws every bar in one colour. `"class"` derives the structural class
#'    with \code{\link{ps_motif_class}}. Alternatively supply a character or
#'    factor vector, either one entry per motif or named by matrix identifier.
#' @param FDR Optional significance cutoff applied before ranking. `NULL`, the
#'    default, keeps every motif.
#'
#' @return A `ggplot` object. Nothing is drawn until the object is printed, so
#'    it can be modified with further `ggplot2` layers first.
#'
#' @details
#' Ranking follows from `statistic` rather than being controlled separately.
#' `"ZSCORE"`, `"FG_AVG"` and `"BG_AVG"` rank from the largest value down.
#' `"P.VALUE"` and `"FDR"` rank from the smallest up, and are plotted as
#' \eqn{-\log_{10}}, since a bar chart of raw p-values spanning many orders of
#' magnitude is unreadable.
#'
#' Bars start at zero, which is meaningful for a z-score and for the
#' \eqn{-\log_{10}} transforms alike.
#'
#' Up to eight groups are drawn from a fixed palette chosen so that adjacent
#' colours remain distinguishable under the common forms of colour-vision
#' deficiency; beyond eight the `ggplot2` default applies. Motif names are
#' always shown on the axis, so the grouping is never carried by colour alone.
#'
#' @seealso \code{\link{ps_motif_class}}, \code{\link{ps_results_table}}
#'
#' @export
#'
#' @importFrom stats setNames
#' @importFrom ggplot2 .data
#'
#' @examples
#' prom_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(prom_path)[1:10]
#' J2020_path <- system.file("extdata", "J2020.rds", package = "PscanR")
#' J2020 <- readRDS(J2020_path)
#' bg_path <- system.file("extdata",
#'     "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'     package = "PscanR"
#' )
#' bg <- ps_retrieve_bg_from_file(bg_path, J2020)
#' bg <- bg[c(
#'     "MA0506.1", "MA0632.2", "MA0611.1",
#'     "MA0685.1", "MA0698.1", "MA0699.1"
#' )]
#' results <- pscan(prom_seq, bg, BPPARAM = BiocParallel::SerialParam())
#'
#' ps_motif_barplot(results, n = 6)
#'
#' # Colour by the structural class of the factor.
#' ps_motif_barplot(results, n = 6, group = "class")
ps_motif_barplot <- function(pfms, n = 20, statistic = c(
                                 "ZSCORE", "P.VALUE", "FDR",
                                 "FG_AVG", "BG_AVG"
                             ),
                             group = NULL, FDR = NULL) {
    statistic <- match.arg(statistic)

    parts <- .ps_barplot_inputs(pfms, group)
    res_table <- parts$table
    grouping <- parts$grouping

    if (!is.null(FDR)) {
    if (!is.numeric(FDR) || length(FDR) != 1L || is.na(FDR)) {
        stop("FDR must be a single numeric value", call. = FALSE)
    }
    keep <- res_table$FDR <= FDR
    res_table <- res_table[keep, , drop = FALSE]
    grouping <- grouping[keep]
    }
    if (nrow(res_table) == 0L) {
    stop("No motifs left to plot after filtering", call. = FALSE)
    }

    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("n must be a single positive number", call. = FALSE)
    }

    # Smaller is better for the two p-value columns, and they are plotted on a
    # log scale; everything else is read directly and ranked downwards.
    ascending <- statistic %in% c("P.VALUE", "FDR")
    value <- res_table[[statistic]]
    ord <- order(value, decreasing = !ascending)
    keep <- utils::head(ord, min(as.integer(n), nrow(res_table)))

    plot_data <- data.frame(
    motif = res_table$NAME[keep],
    motif_id = row.names(res_table)[keep],
    value = if (ascending) -log10(value[keep]) else value[keep],
    stringsAsFactors = FALSE
    )
    # Ties on NAME would collapse bars, so order the factor on the identifier.
    plot_data$label <- factor(
    plot_data$motif_id,
    levels = rev(plot_data$motif_id),
    labels = rev(plot_data$motif)
    )
    axis_label <- if (ascending) {
    paste0("-log10(", statistic, ")")
    } else {
    statistic
    }

    if (is.null(grouping)) {
    plot <- ggplot2::ggplot(
        plot_data, ggplot2::aes(x = .data$value, y = .data$label)
    ) +
        ggplot2::geom_col(width = 0.72, fill = "grey35")
    } else {
    plot_data$group <- grouping[keep]
    n_groups <- length(unique(plot_data$group))
    plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = .data$value, y = .data$label, fill = .data$group)
    ) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::labs(fill = NULL) +
        # Class labels are long, so let the legend wrap rather than clip.
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2))
    if (n_groups <= length(.PS_GROUP_PALETTE)) {
        plot <- plot + ggplot2::scale_fill_manual(
        values = utils::head(.PS_GROUP_PALETTE, n_groups)
        )
    }
    }

    plot +
    ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(x = axis_label, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        legend.justification = "left"
    )
}

# Accepts either a scan result or a results table, and resolves `group` to a
# character vector aligned with the table rows (or NULL).
.ps_barplot_inputs <- function(pfms, group) {
    if (is(pfms, "PFMatrixList")) {
    res_table <- ps_results_table(pfms)
    grouping <- .ps_resolve_group(group, pfms, row.names(res_table))
    return(list(table = res_table, grouping = grouping))
    }
    if (is.data.frame(pfms)) {
    required <- c("NAME", "ZSCORE", "P.VALUE", "FDR")
    missing_cols <- setdiff(required, names(pfms))
    if (length(missing_cols) > 0L) {
        stop(
        "'pfms' looks like a data.frame but is missing column(s): ",
        paste(missing_cols, collapse = ", "),
        ". Pass a PSMatrixList or the output of ps_results_table().",
        call. = FALSE
        )
    }
    if (is.character(group) && length(group) == 1L) {
        stop(
        "group = \"", group, "\" needs motif metadata, which a results ",
        "table does not carry. Pass the PSMatrixList instead, or supply ",
        "'group' as a vector.",
        call. = FALSE
        )
    }
    grouping <- .ps_resolve_group(group, NULL, row.names(pfms))
    return(list(table = pfms, grouping = grouping))
    }
    stop(
    "pfms must be a PSMatrixList or a data.frame from ps_results_table()",
    call. = FALSE
    )
}

.ps_resolve_group <- function(group, pfms, motif_ids) {
    if (is.null(group)) {
    return(NULL)
    }
    if (is.character(group) && length(group) == 1L &&
        group %in% c("class", "family")) {
    if (identical(group, "class")) {
        return(unname(ps_motif_class(pfms)[motif_ids]))
    }
    families <- vapply(pfms, function(x) {
        value <- tags(x)$family
        if (is.null(value) || !nzchar(value[[1L]])) {
        return("Unclassified")
        }
        as.character(value)[[1L]]
    }, character(1L))
    return(unname(stats::setNames(families, ID(pfms))[motif_ids]))
    }
    if (!is.character(group) && !is.factor(group)) {
    stop(
        "group must be NULL, \"class\", \"family\", or a character or ",
        "factor vector",
        call. = FALSE
    )
    }
    group <- as.character(group)
    if (!is.null(names(group))) {
    return(unname(group[motif_ids]))
    }
    if (length(group) != length(motif_ids)) {
    stop(
        "group has ", length(group), " entries but there are ",
        length(motif_ids), " motifs. Supply one entry per motif, or name ",
        "the vector by matrix identifier.",
        call. = FALSE
    )
    }
    group
}

#' Hit Position Against Hit Score
#'
#' Plots every promoter's best hit as a point positioned by where the site sits
#' and by how well it scores, so that position and strength are read together.
#' A positional profile alone cannot show whether a concentration is made of
#' strong sites or weak ones.
#'
#' @param x A `PSMatrix`, or a `PSMatrixList` in which case the motifs are
#'    drawn as panels sharing both axes.
#' @param shift Integer positional shift applied to hit positions, to place
#'    them relative to the TSS. Default `0`.
#'
#'
#' @return A `ggplot` object.
#'
#' @details
#' Two horizontal references are drawn per motif: the background mean, and one
#' standard deviation above it. They are the same thresholds `st` names
#' elsewhere in the package, and they split the points into three bands, each
#' drawn in its own colour: below the mean, between the mean and one standard
#' deviation above it, and beyond that. A point in the top band is a site that
#' `st = "strict"` keeps; the top two bands together are what `"loose"` keeps.
#' Colouring the middle band separately rather than lumping it with the weak
#' hits shows how much of a positional pattern rests on marginal sites.
#'
#' Because each promoter contributes exactly one point, vertical structure is
#' the score distribution and horizontal structure is the positional one. A
#' motif whose strong sites are positionally constrained shows points banking
#' into one region above the upper line while the weak sites stay spread out.
#'
#' @seealso \code{\link{ps_density_plot}}, \code{\link{ps_hitpos_map}}
#'
#' @export
#'
#' @importFrom ggplot2 .data
#'
#' @examples
#' prom_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
#' prom_seq <- readRDS(prom_path)[1:10]
#' J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
#' bg_path <- system.file("extdata",
#'     "J2020_hg38_200u_50d_UCSC.psbg.txt",
#'     package = "PscanR"
#' )
#' bg <- ps_retrieve_bg_from_file(bg_path, J2020)
#' bg <- bg[c("MA0506.1", "MA0632.2")]
#' results <- pscan(prom_seq, bg, BPPARAM = BiocParallel::SerialParam())
#'
#' ps_hit_score_plot(results, shift = -200)
ps_hit_score_plot <- function(x, shift = 0) {
    motifs <- if (is(x, "PFMatrixList")) {
    as.list(x)
    } else if (is(x, "PSMatrix")) {
    list(x)
    } else {
    stop("x must be a PSMatrix or a PSMatrixList", call. = FALSE)
    }

    points <- do.call(rbind, lapply(motifs, function(m) {
    loose <- ps_bg_avg(m)
    strict <- loose + ps_bg_std_dev(m)
    scores <- ps_hits_score(m)
    keep <- !is.na(scores)
    scores <- scores[keep]
    data.frame(
        motif = name(m),
        position = ps_hits_pos(m, pos_shift = shift)[keep],
        score = scores,
        band = .PS_SCORE_BANDS[
        1L + (scores >= loose) + (scores >= strict)
        ],
        stringsAsFactors = FALSE
    )
    }))
    points$band <- factor(points$band, levels = rev(.PS_SCORE_BANDS))
    # Panels follow the order the motifs were given, not the alphabet.
    motif_names <- vapply(motifs, name, character(1L))
    points$motif <- factor(points$motif, levels = unique(motif_names))

    references <- data.frame(
    motif = factor(motif_names, levels = unique(motif_names)),
    loose = vapply(motifs, ps_bg_avg, numeric(1L)),
    strict = vapply(motifs, function(m) {
        ps_bg_avg(m) + ps_bg_std_dev(m)
    }, numeric(1L))
    )

    plot <- ggplot2::ggplot(
    points, ggplot2::aes(x = .data$position, y = .data$score)
    ) +
    ggplot2::geom_hline(
        data = references, ggplot2::aes(yintercept = .data$loose),
        colour = "grey65", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
        data = references, ggplot2::aes(yintercept = .data$strict),
        colour = "grey40", linetype = "dashed", linewidth = 0.4
    ) +
    ggplot2::geom_point(
        ggplot2::aes(colour = .data$band),
        size = 0.9, alpha = 0.55
    ) +
    ggplot2::scale_colour_manual(
        values = stats::setNames(
        c(.PS_DENSITY_COLOUR, .PS_GROUP_PALETTE[[2L]], "grey72"),
        rev(.PS_SCORE_BANDS)
        )
    ) +
    ggplot2::labs(
        x = "Position along promoters", y = "Hit score", colour = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        legend.justification = "left",
        strip.text = ggplot2::element_text(face = "bold")
    )

    if (nlevels(points$motif) > 1L) {
    plot <- plot + ggplot2::facet_wrap(~motif)
    } else {
    plot <- plot + ggplot2::labs(title = levels(points$motif)[[1L]])
    }
    plot
}
