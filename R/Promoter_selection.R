#' Load MANE Select and RefSeq Select Transcript Metadata
#'
#' Downloads and caches the human MANE Select and RefSeq Select transcript
#' tables from UCSC hg38. The resulting table can be passed to
#' \code{\link{ps_select_promoters}} to choose one representative RefSeq
#' promoter per gene.
#'
#' @param organism Character. Currently only \code{"hg38"} is supported.
#' @param cache Logical. If \code{TRUE}, read/write a cached RDS file under
#'   \code{cache_dir}.
#' @param cache_dir Directory used for cached transcript metadata.
#'
#' @return A \code{data.frame} with columns \code{transcript_id},
#'   \code{transcript_base} (the versionless accession), \code{gene_symbol},
#'   \code{selection_source}, and \code{selection_priority}.
#'
#' @details
#' MANE Select is preferred over RefSeq Select. This helper is intentionally
#'   separated from \code{\link{ps_select_promoters}} so callers can cache or
#'   inspect the selection table once and reuse it in multiple analyses.
#'
#' The automatic metadata download provided by this function is human-specific
#'   and applies only to the RefSeq transcript scheme. It currently uses UCSC
#'   hg38 MANE and RefSeq Select tracks. For every other reference, call
#'   \code{\link{ps_select_promoters}} with \code{mode = "representative"},
#'   which applies the ranking rule appropriate to the detected transcript
#'   scheme, or supply your own \code{select_transcripts} table.
#'
#' @seealso \code{\link{ps_select_promoters}}
#'
#' @export
#'
#' @examples
#' cache_dir <- tempfile("pscanr-select-")
#' dir.create(cache_dir)
#' cache_file <- file.path(
#'     cache_dir,
#'     "hg38_mane_refseq_select_transcripts.rds"
#' )
#' example_tx <- data.frame(
#'     transcript_id = c("NM_000546.6", "NM_000492.4"),
#'     transcript_base = c("NM_000546", "NM_000492"),
#'     gene_symbol = c("TP53", "CFTR"),
#'     selection_source = c("MANE Select", "RefSeq Select"),
#'     selection_priority = c(1L, 2L)
#' )
#' saveRDS(example_tx, cache_file)
#'
#' select_tx <- ps_load_select_transcripts(cache_dir = cache_dir)
#' head(select_tx)
#' unlink(cache_dir, recursive = TRUE)
#'
#' if (interactive()) {
#'     current_select_tx <- ps_load_select_transcripts()
#'     head(current_select_tx)
#' }
ps_load_select_transcripts <- function(
    organism = "hg38", cache = TRUE,
    cache_dir = tools::R_user_dir("PscanR", "cache")
) {
    if (!identical(organism, "hg38")) {
        stop("Only organism = 'hg38' is currently supported.")
    }

    cache_file <- file.path(
        cache_dir,
        "hg38_mane_refseq_select_transcripts.rds"
    )
    if (cache && file.exists(cache_file)) {
        return(readRDS(cache_file))
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required to download select transcripts.")
    }

    select_transcripts <- rbind(
        .ps_mane_select_transcripts(),
        .ps_refseq_select_transcripts()
    )
    select_transcripts <- select_transcripts[order(
        select_transcripts$selection_priority,
        select_transcripts$transcript_base
    ), ]
    select_transcripts <- select_transcripts[!duplicated(
        select_transcripts$transcript_base
    ), ]
    row.names(select_transcripts) <- NULL

    if (cache) {
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        saveRDS(select_transcripts, cache_file)
    }
    select_transcripts
}

.ps_mane_select_transcripts <- function() {
    mane <- .ps_ucsc_track_by_chrom("mane")
    mane <- mane[mane$maneStat == "MANE Select" &
        !is.na(mane$ncbiId) & mane$ncbiId != "", ]
    data.frame(
        transcript_id = mane$ncbiId,
        transcript_base = .ps_transcript_base(mane$ncbiId, "refseq"),
        gene_symbol = mane$geneName2,
        selection_source = "MANE Select",
        selection_priority = 1L,
        stringsAsFactors = FALSE
    )
}

.ps_refseq_select_transcripts <- function() {
    refseq <- .ps_ucsc_track_by_chrom("ncbiRefSeqSelect")
    refseq <- refseq[!is.na(refseq$name) & refseq$name != "", ]
    data.frame(
        transcript_id = refseq$name,
        transcript_base = .ps_transcript_base(refseq$name, "refseq"),
        gene_symbol = refseq$name2,
        selection_source = "RefSeq Select",
        selection_priority = 2L,
        stringsAsFactors = FALSE
    )
}

#' Select One Representative Promoter per Gene
#'
#' Maps genes to transcript promoters and selects one promoter per gene using a
#' rule appropriate to the transcript identifier scheme in use. The scheme is
#' detected automatically and reported, so the choice is always auditable.
#'
#' @param genes Character vector of gene symbols, or values in \code{keytype}
#'   if \code{annotation} is not supplied.
#' @param promoter_sequences Optional \code{DNAStringSet} of promoter
#'   sequences, named by transcript identifier.
#' @param promoter_ids Optional character vector of available promoter
#'   transcript IDs. Ignored when \code{promoter_sequences} is supplied.
#' @param annotation Optional \code{data.frame} containing gene-to-transcript
#'   mappings. If omitted, \code{AnnotationDbi::select()} is used with
#'   \code{org_db}.
#' @param org_db Optional organism annotation package object. If omitted and
#'   \code{annotation} is omitted, \code{org.Hs.eg.db::org.Hs.eg.db} is used.
#' @param keytype Key type used for \code{genes} when querying
#'   \code{org_db}. Defaults to \code{"SYMBOL"}.
#' @param gene_col Column in \code{annotation} containing gene identifiers.
#' @param transcript_col Column in \code{annotation} containing transcript
#'   identifiers.
#' @param select_transcripts Optional table from
#'   \code{\link{ps_load_select_transcripts}}. If omitted and
#'   \code{mode = "select"}, the table is downloaded or loaded from cache.
#' @param scheme Transcript identifier scheme. \code{"auto"} (the default)
#'   detects it from the promoter identifiers and reports what it found. Use
#'   \code{"refseq"}, \code{"tair"}, \code{"sgd"} or \code{"generic"} to
#'   override. See Details.
#' @param mode Character. \code{"select"} uses MANE Select, then RefSeq Select,
#'   then the scheme's own ranking; it requires the RefSeq scheme.
#'   \code{"representative"} uses only the scheme's own ranking.
#'   \code{"all_transcripts"} returns every mapped transcript instead of one
#'   per gene.
#' @param fallback Logical. If \code{TRUE}, genes without MANE/RefSeq Select
#'   transcripts fall back to the scheme's own ranking.
#' @param return Character. \code{"mapping"} returns a mapping table,
#'   \code{"sequences"} returns selected promoter sequences, and
#'   \code{"both"} returns a list with both.
#' @param sequence_names Character. How returned sequences should be named:
#'   \code{"gene_transcript"} uses \code{GENE|TRANSCRIPT}; \code{"transcript"}
#'   uses the selected promoter identifier alone.
#' @param quiet Logical. Suppress the scheme-detection message.
#'
#' @return A \code{data.frame}, \code{DNAStringSet}, or list depending on
#'   \code{return}. The mapping table has one row per selected promoter with
#'   columns \code{gene}, \code{transcript_id}, \code{transcript_base},
#'   \code{promoter_id}, \code{scheme}, \code{suffix_role},
#'   \code{transcript_class}, \code{selection_source},
#'   \code{selection_priority}, \code{n_candidates}, \code{decided_by} and
#'   \code{input_order}. Genes that could not be mapped are recorded in
#'   \code{attr(mapping, "ps_unmapped_genes")}; see
#'   \code{\link{ps_selection_summary}}.
#'
#' @details
#' Transcript identifiers do not all mean the same thing, and the difference
#'   matters. In RefSeq, \code{NM_000546.6} carries a release \emph{version}:
#'   two versions name the same transcript, so the suffix is normalised away
#'   before matching. In TAIR, \code{AT1G01110.2} carries a \emph{splice
#'   variant}: the suffix is part of the transcript's identity and is never
#'   stripped, because variants of the same gene frequently have different
#'   transcription start sites. SGD systematic names such as \code{YAL068W-A}
#'   carry no suffix at all.
#'
#' Detected schemes and their ranking rules:
#'   \describe{
#'     \item{\code{"refseq"}}{Human, mouse and Drosophila backgrounds. Prefers
#'       \code{NM_}, then \code{NR_}, then \code{XM_}, then \code{XR_}, then
#'       anything else. Protein accessions (\code{NP_}, \code{XP_}) are not
#'       transcripts and are excluded.}
#'     \item{\code{"tair"}}{Arabidopsis. Prefers the lowest splice-variant
#'       number, which is the primary gene model.}
#'     \item{\code{"sgd"}}{Saccharomyces cerevisiae, which has one transcript
#'       per ORF.}
#'     \item{\code{"generic"}}{Anything unrecognised. Identifiers are matched
#'       whole, nothing is stripped, and all transcripts rank equally. This is
#'       the safe fallback and is chosen with a warning.}
#'   }
#'
#' Ties are broken by a natural ordering of the promoter identifier, in which
#'   embedded numbers compare numerically, so \code{.10} follows \code{.9}.
#'   Because promoter identifiers are unique this is a total order, and the
#'   result therefore does not depend on the order in which promoters were
#'   supplied.
#'
#' Automatic MANE Select and RefSeq Select metadata are available only for
#'   human hg38, so \code{mode = "select"} requires the RefSeq scheme and is an
#'   error otherwise. Use \code{mode = "representative"} elsewhere, or supply a
#'   species-specific \code{select_transcripts} table with
#'   \code{transcript_id} or \code{transcript_base}, \code{selection_source}
#'   and \code{selection_priority} columns.
#'
#' @seealso \code{\link{ps_selection_summary}},
#'   \code{\link{ps_load_select_transcripts}}
#'
#' @export
#'
#' @examples
#' annotation <- data.frame(
#'     SYMBOL = c("GENE1", "GENE1", "GENE2"),
#'     REFSEQ = c("NM_000001.1", "NM_000002.1", "NM_000003.1")
#' )
#' select_tx <- data.frame(
#'     transcript_id = "NM_000002.1",
#'     transcript_base = "NM_000002",
#'     gene_symbol = "GENE1",
#'     selection_source = "MANE Select",
#'     selection_priority = 1L
#' )
#' ps_select_promoters(
#'     c("GENE1", "GENE2"),
#'     annotation = annotation,
#'     select_transcripts = select_tx
#' )
#'
#' # Arabidopsis splice variants are preserved, and the primary model wins.
#' at_annotation <- data.frame(
#'     GENE = c("AT1G01110", "AT1G01110"),
#'     TX = c("AT1G01110.2", "AT1G01110.1")
#' )
#' ps_select_promoters(
#'     "AT1G01110",
#'     promoter_ids = c("AT1G01110.2", "AT1G01110.1"),
#'     annotation = at_annotation,
#'     gene_col = "GENE", transcript_col = "TX",
#'     mode = "representative"
#' )
ps_select_promoters <- function(genes, promoter_sequences = NULL,
                                promoter_ids = NULL, annotation = NULL,
                                org_db = NULL, keytype = "SYMBOL",
                                gene_col = keytype,
                                transcript_col = "REFSEQ",
                                select_transcripts = NULL,
                                scheme = "auto",
                                mode = c(
                                    "select", "representative",
                                    "all_transcripts"
                                ),
                                fallback = TRUE,
                                return = c("mapping", "sequences", "both"),
                                sequence_names = c(
                                    "gene_transcript",
                                    "transcript"
                                ),
                                quiet = FALSE) {
    mode <- match.arg(mode)
    return <- match.arg(return)
    sequence_names <- match.arg(sequence_names)

    input <- .ps_validate_promoter_inputs(
        genes, promoter_sequences, promoter_ids
    )
    genes <- input$genes
    promoter_ids <- input$promoter_ids
    annotation <- .ps_gene_transcript_annotation(
        genes = genes,
        annotation = annotation,
        org_db = org_db,
        keytype = keytype,
        gene_col = gene_col,
        transcript_col = transcript_col
    )

    scheme <- .ps_resolve_scheme(
        scheme = scheme,
        promoter_ids = promoter_ids,
        annotation_ids = annotation$transcript_id,
        quiet = quiet
    )
    annotation$transcript_base <- .ps_transcript_base(
        annotation$transcript_id, scheme
    )

    if (identical(mode, "select")) {
        if (!identical(scheme, "refseq")) {
            stop(
                "mode = \"select\" uses MANE Select and RefSeq Select ",
                "metadata, which exists only for the RefSeq scheme, but the ",
                "identifiers were resolved as \"", scheme, "\". Use ",
                "mode = \"representative\", or supply a species-specific ",
                "'select_transcripts' table.",
                call. = FALSE
            )
        }
        if (is.null(select_transcripts)) {
            select_transcripts <- ps_load_select_transcripts()
        }
    }

    mapping <- .ps_rank_gene_promoters(
        genes = genes,
        annotation = annotation,
        promoter_ids = promoter_ids,
        select_transcripts = select_transcripts,
        scheme = scheme,
        mode = mode,
        fallback = fallback
    )
    attr(mapping, "ps_unmapped_genes") <- setdiff(genes, mapping$gene)
    attr(mapping, "ps_scheme") <- scheme

    .ps_format_promoter_selection(
        mapping, promoter_sequences, return, sequence_names
    )
}

.ps_norm_gene_vector <- function(genes) {
    genes <- trimws(as.character(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    unique(genes)
}

.ps_ucsc_track_by_chrom <- function(track) {
    canonical <- paste0("chr", c(seq_len(22), "X", "Y"))
    rows <- lapply(canonical, function(chrom) {
        url <- sprintf(
            paste0(
                "https://api.genome.ucsc.edu/getData/track?",
                "genome=hg38;track=%s;chrom=%s"
            ),
            track,
            chrom
        )
        response <- jsonlite::fromJSON(url)
        tbl <- response[[track]]
        if (is.null(tbl) || nrow(tbl) == 0L) {
            return(NULL)
        }
        as.data.frame(tbl, stringsAsFactors = FALSE)
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) {
        return(data.frame())
    }
    do.call(rbind, rows)
}

.ps_validate_promoter_inputs <- function(
    genes, promoter_sequences, promoter_ids
) {
    genes <- .ps_norm_gene_vector(genes)
    if (length(genes) == 0L) {
        stop("'genes' must contain at least one non-empty value.")
    }
    if (is.null(promoter_sequences)) {
        return(list(genes = genes, promoter_ids = promoter_ids))
    }
    if (!inherits(promoter_sequences, "DNAStringSet")) {
        stop("'promoter_sequences' must be a DNAStringSet object.")
    }
    promoter_ids <- names(promoter_sequences)
    if (is.null(promoter_ids) || any(promoter_ids == "")) {
        stop("'promoter_sequences' must have transcript identifier names.")
    }
    list(genes = genes, promoter_ids = promoter_ids)
}

.ps_format_promoter_selection <- function(
    mapping, promoter_sequences, return, sequence_names
) {
    if (return == "mapping") {
        return(mapping)
    }
    if (is.null(promoter_sequences)) {
        stop("'promoter_sequences' is required when return is not 'mapping'.")
    }
    selected_sequences <- .ps_subset_promoter_sequences(
        promoter_sequences, mapping, sequence_names
    )
    if (return == "sequences") {
        attr(selected_sequences, "ps_gene_promoter_map") <- mapping
        return(selected_sequences)
    }
    list(mapping = mapping, sequences = selected_sequences)
}

.ps_gene_transcript_annotation <- function(
    genes, annotation, org_db, keytype, gene_col, transcript_col
) {
    if (is.null(annotation)) {
        if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
            stop("Package 'AnnotationDbi' is required when annotation is NULL.")
        }
        if (is.null(org_db)) {
            if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
                stop(
                    "Package 'org.Hs.eg.db' is required when org_db and ",
                    "annotation are NULL."
                )
            }
            org_db <- org.Hs.eg.db::org.Hs.eg.db
        }
        annotation <- AnnotationDbi::select(
            org_db,
            keys = genes,
            keytype = keytype,
            columns = transcript_col
        )
    }

    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
    required_cols <- c(gene_col, transcript_col)
    missing_cols <- setdiff(required_cols, names(annotation))
    if (length(missing_cols) > 0L) {
        stop(
            "'annotation' is missing required column(s): ",
            paste(missing_cols, collapse = ", ")
        )
    }

    annotation <- data.frame(
        gene = as.character(annotation[[gene_col]]),
        transcript_id = as.character(annotation[[transcript_col]]),
        stringsAsFactors = FALSE
    )
    annotation <- annotation[!is.na(annotation$gene) &
        annotation$gene != "" &
        !is.na(annotation$transcript_id) &
        annotation$transcript_id != "", ]
    annotation[!duplicated(annotation), ]
}

.ps_rank_gene_promoters <- function(genes, annotation, promoter_ids,
                                    select_transcripts, scheme, mode,
                                    fallback) {
    annotation <- annotation[annotation$gene %in% genes, ]
    annotation <- .ps_prepare_promoter_candidates(
        annotation, promoter_ids, select_transcripts, scheme, mode
    )
    if (nrow(annotation) == 0L) {
        return(.ps_empty_promoter_mapping())
    }
    annotation <- .ps_assign_promoter_priorities(
        annotation, scheme, mode, fallback
    )
    if (nrow(annotation) == 0L) {
        return(.ps_empty_promoter_mapping())
    }
    .ps_finalize_promoter_mapping(annotation, genes, scheme, mode)
}

.ps_prepare_promoter_candidates <- function(annotation, promoter_ids,
                                            select_transcripts, scheme, mode) {
    # Protein accessions are not transcripts. An OrgDb REFSEQ query returns
    # them alongside mRNAs (NP_034436 next to NM_010306), so drop them before
    # anything else and report how many there were. Counting here rather than
    # after the promoter join gives the number the user can act on.
    n_protein <- sum(
        .ps_transcript_class(annotation$transcript_id, scheme) == "protein"
    )
    if (n_protein > 0L) {
        annotation <- annotation[
            .ps_transcript_class(annotation$transcript_id, scheme) !=
                "protein",
        ]
    }

    if (!is.null(promoter_ids)) {
        promoter_lookup <- .ps_promoter_lookup(promoter_ids, scheme)
        annotation <- merge(annotation, promoter_lookup,
            by = "transcript_base", all = FALSE, sort = FALSE
        )
    } else {
        annotation$promoter_id <- annotation$transcript_id
    }
    annotation$transcript_class <- .ps_transcript_class(
        annotation$promoter_id, scheme
    )

    if (mode == "select") {
        if (is.null(select_transcripts)) {
            select_transcripts <- data.frame()
        }
        select_transcripts <- .ps_norm_select_transcripts(select_transcripts)
        annotation <- merge(annotation, select_transcripts,
            by = "transcript_base", all.x = TRUE, sort = FALSE
        )
    } else {
        annotation$select_gene_symbol <- NA_character_
        annotation$selection_source <- NA_character_
        annotation$selection_priority <- NA_integer_
    }
    attr(annotation, "ps_excluded_protein") <- n_protein
    annotation
}

.ps_assign_promoter_priorities <- function(annotation, scheme, mode,
                                           fallback) {
    excluded_protein <- attr(annotation, "ps_excluded_protein")
    scheme_priority <- .ps_transcript_priority(annotation$promoter_id, scheme)
    scheme_source <- .ps_scheme_source_label(
        annotation$promoter_id, annotation$transcript_class, scheme
    )

    if (mode == "representative" || mode == "all_transcripts") {
        annotation$selection_priority <- scheme_priority
        annotation$selection_source <- scheme_source
    } else {
        missing_select <- is.na(annotation$selection_priority)
        if (fallback) {
            # Select-table hits keep priority 1/2; everything else is ranked by
            # the scheme, offset so it can never outrank a curated selection.
            annotation$selection_priority[missing_select] <-
                scheme_priority[missing_select] + 2L
            annotation$selection_source[missing_select] <-
                scheme_source[missing_select]
        } else {
            annotation <- annotation[!missing_select, ]
        }
    }
    attr(annotation, "ps_excluded_protein") <- excluded_protein
    annotation
}

.ps_scheme_source_label <- function(promoter_id, transcript_class, scheme) {
    if (identical(scheme, "refseq")) {
        return(paste("RefSeq", transcript_class))
    }
    if (identical(scheme, "tair")) {
        return(ifelse(
            transcript_class == "primary",
            "TAIR primary model", "TAIR alternative model"
        ))
    }
    if (identical(scheme, "sgd")) {
        return(rep_len("SGD single transcript", length(promoter_id)))
    }
    rep_len("generic (unranked)", length(promoter_id))
}

.ps_finalize_promoter_mapping <- function(annotation, genes, scheme, mode) {
    excluded_protein <- attr(annotation, "ps_excluded_protein")
    annotation$input_order <- match(annotation$gene, genes)

    # A total order, so nothing is left to the order in which candidates
    # happened to arrive. All keys are integer ranks computed under radix
    # (C-locale) ordering, which keeps the result identical across locales.
    #
    # The first tie-break is scheme-specific. RefSeq keeps its historical
    # lexicographic ordering on the versionless accession, so existing RefSeq
    # analyses are unaffected by this change. Other schemes use the natural
    # ordering, which is what makes TAIR's .10 sort after .9.
    natural_rank <- .ps_dense_rank(.ps_natural_key(annotation$promoter_id))
    tie_rank <- if (identical(scheme, "refseq")) {
        .ps_dense_rank(annotation$transcript_base)
    } else {
        natural_rank
    }
    annotation <- annotation[order(
        annotation$input_order,
        annotation$selection_priority,
        tie_rank,
        natural_rank
    ), ]

    candidates <- table(annotation$gene)
    annotation$n_candidates <- as.integer(candidates[annotation$gene])

    # A gene was decided by tie-break when its best priority is shared.
    best_shared <- vapply(
        split(annotation$selection_priority, annotation$gene),
        function(p) sum(p == min(p)) > 1L, logical(1)
    )
    annotation$decided_by <- ifelse(
        best_shared[annotation$gene], "tie_break", "rule"
    )

    if (mode != "all_transcripts") {
        annotation <- annotation[!duplicated(annotation$gene), ]
    }
    row.names(annotation) <- NULL
    annotation$scheme <- scheme
    annotation$suffix_role <- .ps_scheme_suffix_role(scheme)
    out <- annotation[, c(
        "gene", "transcript_id", "transcript_base", "promoter_id",
        "scheme", "suffix_role", "transcript_class", "selection_source",
        "selection_priority", "n_candidates", "decided_by", "input_order"
    )]
    attr(out, "ps_excluded_protein") <- excluded_protein
    out
}

# Every promoter identifier is retained. Collapsing candidates here, before the
# ranking runs, is what used to make the result depend on input order.
.ps_promoter_lookup <- function(promoter_ids, scheme) {
    promoter_ids <- as.character(promoter_ids)
    data.frame(
        transcript_base = .ps_transcript_base(promoter_ids, scheme),
        promoter_id = promoter_ids,
        stringsAsFactors = FALSE
    )
}

.ps_norm_select_transcripts <- function(select_transcripts) {
    if (nrow(select_transcripts) == 0L) {
        return(data.frame(
            transcript_base = character(),
            select_gene_symbol = character(),
            selection_source = character(),
            selection_priority = integer()
        ))
    }
    select_transcripts <- as.data.frame(select_transcripts,
        stringsAsFactors = FALSE
    )
    if (!"transcript_base" %in% names(select_transcripts)) {
        if (!"transcript_id" %in% names(select_transcripts)) {
            stop(
                "'select_transcripts' must contain 'transcript_base' or ",
                "'transcript_id'."
            )
        }
        select_transcripts$transcript_base <- .ps_transcript_base(
            select_transcripts$transcript_id, "refseq"
        )
    }
    if (!"selection_source" %in% names(select_transcripts)) {
        select_transcripts$selection_source <- "Select transcript"
    }
    if (!"selection_priority" %in% names(select_transcripts)) {
        select_transcripts$selection_priority <- 1L
    }
    if (!"gene_symbol" %in% names(select_transcripts)) {
        select_transcripts$gene_symbol <- NA_character_
    }

    out <- data.frame(
        transcript_base = as.character(select_transcripts$transcript_base),
        select_gene_symbol = as.character(select_transcripts$gene_symbol),
        selection_source = as.character(select_transcripts$selection_source),
        selection_priority = as.integer(select_transcripts$selection_priority),
        stringsAsFactors = FALSE
    )
    out <- out[order(out$selection_priority, out$transcript_base), ]
    out[!duplicated(out$transcript_base), ]
}

.ps_empty_promoter_mapping <- function() {
    data.frame(
        gene = character(),
        transcript_id = character(),
        transcript_base = character(),
        promoter_id = character(),
        scheme = character(),
        suffix_role = character(),
        transcript_class = character(),
        selection_source = character(),
        selection_priority = integer(),
        n_candidates = integer(),
        decided_by = character(),
        input_order = integer()
    )
}

.ps_subset_promoter_sequences <- function(
    promoter_sequences, mapping, sequence_names
) {
    seqs <- promoter_sequences[mapping$promoter_id]
    if (sequence_names == "gene_transcript") {
        names(seqs) <- paste(mapping$gene, mapping$promoter_id, sep = "|")
    } else {
        names(seqs) <- mapping$promoter_id
    }
    seqs
}

#' Summarise a Promoter Selection
#'
#' Reports how a mapping produced by \code{\link{ps_select_promoters}} was
#' arrived at: which scheme was used, which rule selected each promoter, how
#' many genes had a real choice to make, and how many were decided by the
#' tie-break rather than by the ranking.
#'
#' @param mapping A mapping \code{data.frame} returned by
#'   \code{\link{ps_select_promoters}}, or the list returned with
#'   \code{return = "both"}.
#'
#' @return A list with two elements. \code{overall} is a one-row
#'   \code{data.frame} of counts: genes requested, genes mapped, genes
#'   unmapped, genes with more than one candidate, genes decided by tie-break,
#'   and candidates excluded as protein accessions. \code{by_source} is a
#'   \code{data.frame} with one row per \code{selection_source}.
#'
#' @seealso \code{\link{ps_select_promoters}}
#'
#' @export
#'
#' @examples
#' annotation <- data.frame(
#'     GENE = c("AT1G01110", "AT1G01110", "AT1G01160"),
#'     TX = c("AT1G01110.2", "AT1G01110.1", "AT1G01160.1")
#' )
#' mapping <- ps_select_promoters(
#'     c("AT1G01110", "AT1G01160"),
#'     promoter_ids = c("AT1G01110.2", "AT1G01110.1", "AT1G01160.1"),
#'     annotation = annotation,
#'     gene_col = "GENE", transcript_col = "TX",
#'     mode = "representative"
#' )
#' ps_selection_summary(mapping)
ps_selection_summary <- function(mapping) {
    if (is.list(mapping) && !is.data.frame(mapping) &&
        !is.null(mapping$mapping)) {
        mapping <- mapping$mapping
    }
    if (!is.data.frame(mapping) || !"promoter_id" %in% names(mapping)) {
        stop(
            "'mapping' must be a promoter selection returned by ",
            "ps_select_promoters().",
            call. = FALSE
        )
    }
    unmapped <- attr(mapping, "ps_unmapped_genes")
    if (is.null(unmapped)) {
        unmapped <- character()
    }
    excluded <- attr(mapping, "ps_excluded_protein")
    if (is.null(excluded)) {
        excluded <- NA_integer_
    }
    overall <- data.frame(
        scheme = if (nrow(mapping) > 0L) {
            mapping$scheme[[1L]]
        } else {
            NA_character_
        },
        genes_requested = length(unique(mapping$gene)) + length(unmapped),
        genes_mapped = length(unique(mapping$gene)),
        genes_unmapped = length(unmapped),
        genes_with_choice = length(unique(
            mapping$gene[mapping$n_candidates > 1L]
        )),
        genes_decided_by_tie_break = length(unique(
            mapping$gene[mapping$decided_by == "tie_break"]
        )),
        candidates_excluded_protein = as.integer(excluded),
        stringsAsFactors = FALSE
    )
    counts <- table(mapping$selection_source)
    by_source <- data.frame(
        selection_source = names(counts),
        promoters = as.integer(counts),
        stringsAsFactors = FALSE
    )
    by_source <- by_source[order(-by_source$promoters), ]
    row.names(by_source) <- NULL
    list(overall = overall, by_source = by_source)
}
