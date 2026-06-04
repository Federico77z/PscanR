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
#' @return A \code{data.frame} with RefSeq transcript IDs, versionless RefSeq
#'   IDs, gene symbols, selection source, and selection priority.
#'
#' @details
#' MANE Select is preferred over RefSeq Select. This helper is intentionally
#'   separated from \code{\link{ps_select_promoters}} so callers can cache or
#'   inspect the selection table once and reuse it in multiple analyses.
#'
#' The automatic metadata download provided by this function is human-specific.
#'   It currently uses UCSC hg38 MANE and RefSeq Select tracks. For mouse,
#'   Drosophila melanogaster, Saccharomyces cerevisiae, or other organisms,
#'   users should either call \code{\link{ps_select_promoters}} with
#'   \code{mode = "representative"} or provide a species-specific
#'   \code{select_transcripts} table to \code{\link{ps_select_promoters}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' select_tx <- ps_load_select_transcripts()
#' head(select_tx)
#' }
ps_load_select_transcripts <- function(organism = "hg38", cache = TRUE,
                                       cache_dir = tools::R_user_dir(
                                        "PscanR", "cache"
                                       )) {
    if (!identical(organism, "hg38")) {
        stop("Only organism = 'hg38' is currently supported.")
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required to download select transcripts.")
    }

    cache_file <- file.path(
        cache_dir,
        "hg38_mane_refseq_select_transcripts.rds"
    )
    if (cache && file.exists(cache_file)) {
        return(readRDS(cache_file))
    }

    mane <- .ps_ucsc_track_by_chrom("mane")
    mane <- mane[mane$maneStat == "MANE Select" &
        !is.na(mane$ncbiId) & mane$ncbiId != "", ]
    mane_select <- data.frame(
        refseq_id = mane$ncbiId,
        refseq_clean = .ps_clean_refseq_id(mane$ncbiId),
        gene_symbol = mane$geneName2,
        selection_source = "MANE Select",
        selection_priority = 1L,
        stringsAsFactors = FALSE
    )

    refseq <- .ps_ucsc_track_by_chrom("ncbiRefSeqSelect")
    refseq <- refseq[!is.na(refseq$name) & refseq$name != "", ]
    refseq_select <- data.frame(
        refseq_id = refseq$name,
        refseq_clean = .ps_clean_refseq_id(refseq$name),
        gene_symbol = refseq$name2,
        selection_source = "RefSeq Select",
        selection_priority = 2L,
        stringsAsFactors = FALSE
    )

    select_transcripts <- rbind(mane_select, refseq_select)
    select_transcripts <- select_transcripts[order(
        select_transcripts$selection_priority,
        select_transcripts$refseq_clean
    ), ]
    select_transcripts <- select_transcripts[!duplicated(
        select_transcripts$refseq_clean
    ), ]
    row.names(select_transcripts) <- NULL

    if (cache) {
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        saveRDS(select_transcripts, cache_file)
    }
    select_transcripts
}

#' Select One Representative Promoter per Gene
#'
#' Maps gene symbols to RefSeq transcript promoters and selects one promoter
#' per gene by preferring MANE Select, then RefSeq Select, then a deterministic
#' RefSeq fallback.
#'
#' @param genes Character vector of gene symbols, or values in \code{keytype}
#'   if \code{annotation} is not supplied.
#' @param promoter_sequences Optional \code{DNAStringSet} of promoter
#'   sequences. Names must be RefSeq transcript IDs, with or without versions.
#' @param promoter_ids Optional character vector of available promoter
#'   transcript IDs. Ignored when \code{promoter_sequences} is supplied.
#' @param annotation Optional \code{data.frame} containing gene-to-RefSeq
#'   transcript mappings. If omitted, \code{AnnotationDbi::select()} is used
#'   with \code{org_db}.
#' @param org_db Optional organism annotation package object. If omitted and
#'   \code{annotation} is omitted, \code{org.Hs.eg.db::org.Hs.eg.db} is used.
#' @param keytype Key type used for \code{genes} when querying
#'   \code{org_db}. Defaults to \code{"SYMBOL"}.
#' @param gene_col Column in \code{annotation} containing gene identifiers.
#' @param refseq_col Column in \code{annotation} containing RefSeq IDs.
#' @param select_transcripts Optional table from
#'   \code{\link{ps_load_select_transcripts}}. If omitted and
#'   \code{mode = "select"}, the table is downloaded or loaded from cache.
#' @param mode Character. \code{"select"} uses MANE Select, RefSeq Select,
#'   then RefSeq fallback. \code{"representative"} ignores MANE/RefSeq Select
#'   and uses the deterministic RefSeq fallback. \code{"all_transcripts"}
#'   returns every mapped transcript instead of one per gene.
#' @param fallback Logical. If \code{TRUE}, genes without MANE/RefSeq Select
#'   transcripts can use the deterministic RefSeq fallback.
#' @param return Character. \code{"mapping"} returns a mapping table,
#'   \code{"sequences"} returns selected promoter sequences, and
#'   \code{"both"} returns a list with both.
#' @param sequence_names Character. How returned sequences should be named:
#'   \code{"gene_refseq"} uses \code{GENE|REFSEQ}; \code{"refseq"} uses the
#'   selected transcript ID.
#'
#' @return A \code{data.frame}, \code{DNAStringSet}, or list depending on
#'   \code{return}. The mapping table includes the selected RefSeq ID, promoter
#'   ID, selection source, and selection priority.
#'
#' @details
#' The deterministic fallback prefers \code{NM_} transcripts, then
#'   \code{NR_} transcripts, then any other RefSeq IDs; ties are resolved by
#'   lexicographic transcript ID. This mirrors PscanR workflows where promoter
#'   sequence names are RefSeq transcript IDs and keeps one promoter per gene
#'   unless \code{mode = "all_transcripts"} is requested.
#'
#' Automatic MANE Select and RefSeq Select metadata are currently available
#'   only for human hg38 through \code{\link{ps_load_select_transcripts}}.
#'   For non-human organisms supported by PscanR, such as mouse,
#'   Drosophila melanogaster, and Saccharomyces cerevisiae, the function can
#'   still be used in two ways. First, users can set
#'   \code{mode = "representative"} and provide organism-appropriate
#'   \code{annotation} and promoter IDs or sequences; this uses only the
#'   deterministic transcript-ID fallback and is most meaningful for RefSeq-like
#'   transcript IDs. Second, users can provide their own species-specific
#'   \code{select_transcripts} table with \code{refseq_id} or
#'   \code{refseq_clean}, \code{selection_source}, and
#'   \code{selection_priority} columns. This is the recommended route for
#'   Ensembl, FlyBase, SGD, or other non-RefSeq transcript naming schemes.
#'
#' @export
#'
#' @examples
#' annotation <- data.frame(
#'     SYMBOL = c("GENE1", "GENE1", "GENE2"),
#'     REFSEQ = c("NM_000001.1", "NM_000002.1", "NM_000003.1")
#' )
#' select_tx <- data.frame(
#'     refseq_id = "NM_000002.1",
#'     refseq_clean = "NM_000002",
#'     gene_symbol = "GENE1",
#'     selection_source = "MANE Select",
#'     selection_priority = 1L
#' )
#' ps_select_promoters(
#'     c("GENE1", "GENE2"),
#'     annotation = annotation,
#'     select_transcripts = select_tx
#' )
ps_select_promoters <- function(genes, promoter_sequences = NULL,
                                promoter_ids = NULL, annotation = NULL,
                                org_db = NULL, keytype = "SYMBOL",
                                gene_col = keytype,
                                refseq_col = "REFSEQ",
                                select_transcripts = NULL,
                                mode = c("select", "representative",
                                         "all_transcripts"),
                                fallback = TRUE,
                                return = c("mapping", "sequences", "both"),
                                sequence_names = c("gene_refseq",
                                                   "refseq")) {
    mode <- match.arg(mode)
    return <- match.arg(return)
    sequence_names <- match.arg(sequence_names)

    genes <- .ps_norm_gene_vector(genes)
    if (length(genes) == 0L) {
        stop("'genes' must contain at least one non-empty value.")
    }

    if (!is.null(promoter_sequences)) {
        if (!inherits(promoter_sequences, "DNAStringSet")) {
            stop("'promoter_sequences' must be a DNAStringSet object.")
        }
        promoter_ids <- names(promoter_sequences)
        if (is.null(promoter_ids) || any(promoter_ids == "")) {
            stop("'promoter_sequences' must have RefSeq transcript names.")
        }
    }

    annotation <- .ps_gene_refseq_annotation(
        genes = genes,
        annotation = annotation,
        org_db = org_db,
        keytype = keytype,
        gene_col = gene_col,
        refseq_col = refseq_col
    )

    if (mode == "select" && is.null(select_transcripts)) {
        select_transcripts <- ps_load_select_transcripts()
    }

    mapping <- .ps_rank_gene_promoters(
        genes = genes,
        annotation = annotation,
        promoter_ids = promoter_ids,
        select_transcripts = select_transcripts,
        mode = mode,
        fallback = fallback
    )

    if (return == "mapping") {
        return(mapping)
    }
    if (is.null(promoter_sequences)) {
        stop("'promoter_sequences' is required when return is not 'mapping'.")
    }

    selected_sequences <- .ps_subset_promoter_sequences(
        promoter_sequences,
        mapping,
        sequence_names
    )
    if (return == "sequences") {
        attr(selected_sequences, "ps_gene_promoter_map") <- mapping
        return(selected_sequences)
    }
    list(mapping = mapping, sequences = selected_sequences)
}

.ps_clean_refseq_id <- function(x) {
    sub("\\..*$", "", as.character(x))
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

.ps_gene_refseq_annotation <- function(genes, annotation, org_db, keytype,
                                       gene_col, refseq_col) {
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
            columns = refseq_col
        )
    }

    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
    required_cols <- c(gene_col, refseq_col)
    missing_cols <- setdiff(required_cols, names(annotation))
    if (length(missing_cols) > 0L) {
        stop(
            "'annotation' is missing required column(s): ",
            paste(missing_cols, collapse = ", ")
        )
    }

    annotation <- data.frame(
        gene = as.character(annotation[[gene_col]]),
        refseq_id = as.character(annotation[[refseq_col]]),
        stringsAsFactors = FALSE
    )
    annotation <- annotation[!is.na(annotation$gene) &
        annotation$gene != "" &
        !is.na(annotation$refseq_id) &
        annotation$refseq_id != "", ]
    annotation$refseq_clean <- .ps_clean_refseq_id(annotation$refseq_id)
    annotation <- annotation[!duplicated(annotation), ]
    annotation
}

.ps_rank_gene_promoters <- function(genes, annotation, promoter_ids,
                                    select_transcripts, mode, fallback) {
    annotation <- annotation[annotation$gene %in% genes, ]
    if (!is.null(promoter_ids)) {
        promoter_lookup <- .ps_promoter_lookup(promoter_ids)
        annotation <- merge(annotation, promoter_lookup,
            by = "refseq_clean", all = FALSE, sort = FALSE
        )
    } else {
        annotation$promoter_id <- annotation$refseq_clean
    }

    if (nrow(annotation) == 0L) {
        return(.ps_empty_promoter_mapping())
    }

    if (mode == "select") {
        if (is.null(select_transcripts)) {
            select_transcripts <- data.frame()
        }
        select_transcripts <- .ps_norm_select_transcripts(select_transcripts)
        annotation <- merge(annotation, select_transcripts,
            by = "refseq_clean", all.x = TRUE, sort = FALSE
        )
    } else {
        annotation$select_gene_symbol <- NA_character_
        annotation$selection_source <- NA_character_
        annotation$selection_priority <- NA_integer_
    }

    annotation$refseq_type <- .ps_refseq_type(annotation$refseq_clean)
    fallback_priority <- .ps_fallback_priority(annotation$refseq_clean)
    fallback_source <- paste("RefSeq", annotation$refseq_type, "fallback")
    fallback_source[annotation$refseq_type == "other"] <- "RefSeq fallback"

    if (mode == "representative") {
        annotation$selection_priority <- fallback_priority - 2L
        annotation$selection_source <- fallback_source
    } else if (mode == "select") {
        missing_select <- is.na(annotation$selection_priority)
        if (fallback) {
            annotation$selection_priority[missing_select] <-
                fallback_priority[missing_select]
            annotation$selection_source[missing_select] <-
                fallback_source[missing_select]
        } else {
            annotation <- annotation[!missing_select, ]
        }
    } else {
        annotation$selection_priority <- fallback_priority
        annotation$selection_source[is.na(annotation$selection_source)] <-
            fallback_source[is.na(annotation$selection_source)]
    }

    if (nrow(annotation) == 0L) {
        return(.ps_empty_promoter_mapping())
    }

    annotation$input_order <- match(annotation$gene, genes)
    annotation <- annotation[order(
        annotation$input_order,
        annotation$selection_priority,
        annotation$refseq_clean
    ), ]

    if (mode != "all_transcripts") {
        annotation <- annotation[!duplicated(annotation$gene), ]
    }
    row.names(annotation) <- NULL
    annotation[, c(
        "gene", "refseq_id", "refseq_clean", "promoter_id",
        "selection_source", "selection_priority", "refseq_type",
        "input_order"
    )]
}

.ps_promoter_lookup <- function(promoter_ids) {
    promoter_ids <- as.character(promoter_ids)
    data.frame(
        refseq_clean = .ps_clean_refseq_id(promoter_ids),
        promoter_id = promoter_ids,
        stringsAsFactors = FALSE
    )[!duplicated(.ps_clean_refseq_id(promoter_ids)), ]
}

.ps_norm_select_transcripts <- function(select_transcripts) {
    if (nrow(select_transcripts) == 0L) {
        return(data.frame(
            refseq_clean = character(),
            select_gene_symbol = character(),
            selection_source = character(),
            selection_priority = integer()
        ))
    }
    select_transcripts <- as.data.frame(select_transcripts,
        stringsAsFactors = FALSE
    )
    if (!"refseq_clean" %in% names(select_transcripts)) {
        if (!"refseq_id" %in% names(select_transcripts)) {
            stop(
                "'select_transcripts' must contain 'refseq_clean' or ",
                "'refseq_id'."
            )
        }
        select_transcripts$refseq_clean <- .ps_clean_refseq_id(
            select_transcripts$refseq_id
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
        refseq_clean = as.character(select_transcripts$refseq_clean),
        select_gene_symbol = as.character(select_transcripts$gene_symbol),
        selection_source = as.character(select_transcripts$selection_source),
        selection_priority = as.integer(select_transcripts$selection_priority),
        stringsAsFactors = FALSE
    )
    out <- out[order(out$selection_priority, out$refseq_clean), ]
    out[!duplicated(out$refseq_clean), ]
}

.ps_refseq_type <- function(refseq_clean) {
    ifelse(grepl("^NM_", refseq_clean), "NM",
        ifelse(grepl("^NR_", refseq_clean), "NR", "other")
    )
}

.ps_fallback_priority <- function(refseq_clean) {
    ifelse(grepl("^NM_", refseq_clean), 3L,
        ifelse(grepl("^NR_", refseq_clean), 4L, 5L)
    )
}

.ps_empty_promoter_mapping <- function() {
    data.frame(
        gene = character(),
        refseq_id = character(),
        refseq_clean = character(),
        promoter_id = character(),
        selection_source = character(),
        selection_priority = integer(),
        refseq_type = character(),
        input_order = integer()
    )
}

.ps_subset_promoter_sequences <- function(promoter_sequences, mapping,
                                          sequence_names) {
    seqs <- promoter_sequences[mapping$promoter_id]
    if (sequence_names == "gene_refseq") {
        names(seqs) <- paste(mapping$gene, mapping$refseq_clean, sep = "|")
    } else {
        names(seqs) <- mapping$refseq_clean
    }
    seqs
}
