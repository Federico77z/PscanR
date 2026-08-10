# Transcript identifier schemes.
#
# PscanR matches user-supplied transcript identifiers against the promoter
# identifiers carried by a background. Doing that correctly requires knowing
# what the identifier grammar means, because the same punctuation carries
# different semantics in different references:
#
#   NM_000546.6   RefSeq: the suffix is a VERSION. Two versions name the same
#                 transcript, so it must be normalised away before matching.
#   AT1G01110.2   TAIR:   the suffix is a SPLICE VARIANT. Stripping it merges
#                 distinct transcripts onto their shared gene, and for
#                 Arabidopsis most multi-variant genes have variants at
#                 different transcription start sites.
#   YAL068W-A     SGD:    there is no suffix at all; the dash belongs to the
#                 systematic name.
#
# Every identifier normalisation in the package goes through this file. No
# other code should strip an identifier itself.

# Each scheme describes one grammar:
#   pattern      recogniser used for detection
#   base         map an identifier to its stable TRANSCRIPT key -- never a gene
#                key. Only "refseq" removes anything.
#   suffix_role  what the trailing component means: version, variant or none
#   label        human-readable name used in messages
.PS_TRANSCRIPT_SCHEMES <- c("refseq", "sgd", "tair", "generic")

.ps_scheme_pattern <- function(scheme) {
    switch(scheme,
        refseq = "^[A-Z]{2}_[0-9]+(\\.[0-9]+)?$",
        sgd = "^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$",
        tair = "^AT[0-9CM]G[0-9]+\\.[0-9]+$",
        generic = ".",
        stop("Unknown transcript scheme: ", scheme, call. = FALSE)
    )
}

.ps_scheme_label <- function(scheme) {
    switch(scheme,
        refseq = "RefSeq accession with version suffix",
        sgd = "SGD systematic name",
        tair = "TAIR AGI locus with splice-variant suffix",
        generic = "unrecognised identifier grammar",
        scheme
    )
}

.ps_scheme_suffix_role <- function(scheme) {
    switch(scheme,
        refseq = "version",
        sgd = "none",
        tair = "variant",
        generic = "none",
        "none"
    )
}

#' Reduce transcript identifiers to their stable transcript key
#'
#' Only the RefSeq scheme removes anything: its trailing \code{.N} is a release
#' version of the same transcript. For every other scheme the identifier is
#' already a transcript key and is returned unchanged, so distinct splice
#' variants are never merged.
#'
#' @param x Character vector of transcript identifiers.
#' @param scheme One of \code{"refseq"}, \code{"sgd"}, \code{"tair"},
#'   \code{"generic"}.
#'
#' @return Character vector of transcript keys.
#'
#' @noRd
.ps_transcript_base <- function(x, scheme) {
    x <- as.character(x)
    if (identical(scheme, "refseq")) {
        return(sub("\\.[0-9]+$", "", x))
    }
    x
}

#' Classify a transcript identifier within its scheme
#'
#' @return For RefSeq, one of \code{"NM"}, \code{"NR"}, \code{"XM"},
#'   \code{"XR"}, \code{"protein"} or \code{"other"}. For TAIR,
#'   \code{"primary"} for the \code{.1} gene model and \code{"alternative"}
#'   otherwise. For SGD, \code{"single"}.
#'
#' @noRd
.ps_transcript_class <- function(x, scheme) {
    x <- as.character(x)
    if (identical(scheme, "refseq")) {
        return(ifelse(grepl("^NM_", x), "NM",
            ifelse(grepl("^NR_", x), "NR",
                ifelse(grepl("^XM_", x), "XM",
                    ifelse(grepl("^XR_", x), "XR",
                        ifelse(grepl("^[NXY]P_", x), "protein", "other")
                    )
                )
            )
        ))
    }
    if (identical(scheme, "tair")) {
        variant <- .ps_tair_variant(x)
        return(ifelse(is.na(variant), "other",
            ifelse(variant == 1L, "primary", "alternative")
        ))
    }
    if (identical(scheme, "sgd")) {
        return(rep_len("single", length(x)))
    }
    rep_len("other", length(x))
}

.ps_tair_variant <- function(x) {
    suffix <- sub("^.*\\.", "", as.character(x))
    suppressWarnings(as.integer(suffix))
}

#' Rank candidate transcripts within a gene
#'
#' Lower is better. RefSeq prefers curated mRNA over curated non-coding over
#' predicted models; TAIR prefers the lowest splice-variant number, which is
#' the primary gene model; SGD has a single transcript per ORF.
#'
#' Protein accessions get \code{NA}: they are not transcripts and are excluded
#' by the caller rather than ranked.
#'
#' @noRd
.ps_transcript_priority <- function(x, scheme) {
    if (identical(scheme, "refseq")) {
        cls <- .ps_transcript_class(x, scheme)
        return(unname(c(
            NM = 1L, NR = 2L, XM = 3L, XR = 4L,
            protein = NA_integer_, other = 5L
        )[cls]))
    }
    if (identical(scheme, "tair")) {
        variant <- .ps_tair_variant(x)
        # Splice-variant number is itself the rank: .1 is the primary model.
        return(ifelse(is.na(variant), 9999L, variant))
    }
    rep_len(1L, length(x))
}

#' Order identifiers so that numeric components compare numerically
#'
#' Plain lexicographic ordering places \code{AT4G32850.10} before
#' \code{AT4G32850.2}. This zero-pads every embedded digit run to a fixed width
#' so that a byte comparison agrees with a numeric one, then sorts with the
#' radix method, which ignores the locale. The result is a deterministic total
#' order over unique identifiers, identical on every platform.
#'
#' @noRd
.PS_NATURAL_ORDER_WIDTH <- 12L

#' Sort key in which embedded digit runs compare numerically
#'
#' @noRd
.ps_natural_key <- function(x) {
    x <- as.character(x)
    # Delimit digit runs, then widen each one to a fixed number of digits. A
    # run of length L reaches the target width after WIDTH - L passes, so
    # WIDTH - 1 passes always suffice.
    pattern <- sprintf("<([0-9]{1,%d})>", .PS_NATURAL_ORDER_WIDTH - 1L)
    padded <- gsub("([0-9]+)", "<\\1>", x)
    for (i in seq_len(.PS_NATURAL_ORDER_WIDTH - 1L)) {
        padded <- gsub(pattern, "<0\\1>", padded)
    }
    padded
}

.ps_natural_order <- function(x) {
    order(.ps_natural_key(x), method = "radix")
}

#' Dense rank under C-locale byte ordering
#'
#' Equal values receive equal ranks, so a rank never manufactures a
#' distinction that the underlying value does not make. \code{method = "radix"}
#' keeps the ordering independent of the session locale.
#'
#' @noRd
.ps_dense_rank <- function(x) {
    match(x, sort(unique(as.character(x)), method = "radix"))
}

#' Detect the transcript identifier scheme of a character vector
#'
#' @param x Character vector of identifiers.
#'
#' @return A list with \code{scheme}, \code{matched} (number of identifiers
#'   matching the winning pattern), \code{total}, and \code{ambiguous} (a
#'   character vector of competing schemes when detection was not clean).
#'
#' @noRd
.ps_detect_scheme <- function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & x != ""]
    total <- length(x)
    if (total == 0L) {
        return(list(
            scheme = "generic", matched = 0L, total = 0L,
            ambiguous = character()
        ))
    }
    named <- setdiff(.PS_TRANSCRIPT_SCHEMES, "generic")
    hits <- vapply(named, function(scheme) {
        sum(grepl(.ps_scheme_pattern(scheme), x))
    }, integer(1))
    # A scheme wins only if it accounts for every identifier. Anything less is
    # a mixed or unknown set, which must fall back to the lossless grammar.
    complete <- names(hits)[hits == total]
    if (length(complete) == 1L) {
        return(list(
            scheme = complete, matched = total, total = total,
            ambiguous = character()
        ))
    }
    list(
        scheme = "generic",
        matched = if (length(hits)) max(hits) else 0L,
        total = total,
        ambiguous = if (length(complete) > 1L) {
            complete
        } else {
            names(hits)[hits > 0L]
        }
    )
}

#' Resolve the scheme to use, reporting what was decided and why
#'
#' @param scheme Requested scheme, or \code{"auto"}.
#' @param promoter_ids Identifiers from the promoter set. Authoritative.
#' @param annotation_ids Identifiers supplied by the user, checked for
#'   agreement with the promoter set.
#' @param quiet Suppress the informational message.
#'
#' @noRd
.ps_resolve_scheme <- function(scheme = "auto", promoter_ids = NULL,
                               annotation_ids = NULL, quiet = FALSE) {
    scheme <- match.arg(
        scheme, c("auto", .PS_TRANSCRIPT_SCHEMES)
    )
    reference <- if (length(promoter_ids) > 0L) {
        promoter_ids
    } else {
        annotation_ids
    }
    detected <- .ps_detect_scheme(reference)

    if (!identical(scheme, "auto")) {
        if (!quiet) {
            message(sprintf(
                "Using transcript scheme \"%s\" (%s); suffix treated as %s.",
                scheme, .ps_scheme_label(scheme),
                .ps_scheme_suffix_role(scheme)
            ))
        }
        return(scheme)
    }

    resolved <- detected$scheme
    if (identical(resolved, "generic")) {
        warning(
            "Could not identify the transcript identifier scheme",
            if (length(detected$ambiguous) > 0L) {
                paste0(
                    " (partial matches: ",
                    paste(detected$ambiguous, collapse = ", "), ")"
                )
            } else {
                ""
            },
            ". Falling back to scheme \"generic\": identifiers are matched ",
            "whole, nothing is stripped, and every transcript ranks equally. ",
            "Pass 'scheme' explicitly if this is wrong.",
            call. = FALSE
        )
    } else if (!quiet) {
        candidates <- reference[!is.na(reference) &
            grepl(.ps_scheme_pattern(resolved), reference)]
        example <- if (length(candidates) > 0L) candidates[[1L]] else NA_character_
        role <- .ps_scheme_suffix_role(resolved)
        message(sprintf(
            "Detected transcript scheme: \"%s\" (%d/%d identifiers matched)\n  %s -> transcript %s\n  %s",
            resolved, detected$matched, detected$total,
            example, .ps_transcript_base(example, resolved),
            switch(role,
                version = paste(
                    "the trailing .N is a VERSION and is normalised away",
                    "for matching"
                ),
                variant = paste(
                    "the trailing .N is a SPLICE VARIANT and is preserved,",
                    "not stripped"
                ),
                "these identifiers carry no version or variant suffix"
            )
        ))
    }

    # The user's identifiers must speak the same grammar as the promoter set,
    # otherwise nothing will join and the reason would be invisible.
    if (length(promoter_ids) > 0L && length(annotation_ids) > 0L &&
        !identical(resolved, "generic")) {
        annotation_detected <- .ps_detect_scheme(annotation_ids)
        if (!annotation_detected$scheme %in% c(resolved, "generic")) {
            stop(
                "Transcript identifier schemes disagree: the promoter set ",
                "looks like \"", resolved, "\" but 'annotation' looks like \"",
                annotation_detected$scheme, "\". Supply matching identifiers ",
                "or set 'scheme' explicitly.",
                call. = FALSE
            )
        }
    }
    resolved
}
