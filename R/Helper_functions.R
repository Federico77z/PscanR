#' @keywords internal
#' @importFrom TFBSTools PFMatrixList
#'
.ps_warn_missing_bg <- function(pfms) {
    nabg <- vapply(pfms, ps_bg_avg, FUN.VALUE = numeric(length = 1L))
    nabg2 <- vapply(pfms, ps_bg_std_dev, FUN.VALUE = numeric(length = 1L))
    nabg <- is.na(nabg) | is.na(nabg2)
    if (any(nabg)) {
    warning(sprintf(
        "\nNo Pscan background for %s %s", name(pfms)[nabg],
        ID(pfms)[nabg]
    ))
    }
}

#' Warn when the foreground is a large share of the background
#'
#' Why this matters, and the option that moves the threshold, are documented in
#' the "Foreground size" section of \code{\link{pscan}}. That is the copy a
#' user can read, so keep the two in step; the warning names it.
#'
#' @param n Number of foreground sequences.
#' @param bg_sizes Integer vector of per-motif background sizes.
#' @param when Either `"scan"` or `"results"`, selecting the opening clause so
#'    that the two call sites do not read as the same sentence twice.
#'
#' @return `NULL`, invisibly. Called for the warning.
#'
#' @noRd
#'
.ps_warn_foreground_fraction <- function(n, bg_sizes, when = "scan") {
    threshold <- getOption("PscanR.foreground.max_fraction", 0.1)
    # A threshold of zero or less is not a stricter check, it is an
    # unconditional warning: `fraction <= threshold` is then false for every
    # real input. `Inf`, which silences the check, is deliberately allowed.
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold <= 0) {
    threshold <- 0.1
    }
    # Per-motif and NA for any motif absent from the background table, so the
    # median of the usable values stands in for the run.
    sizes <- bg_sizes[!is.na(bg_sizes) & bg_sizes > 0]
    if (length(n) != 1L || is.na(n) || n <= 0 || length(sizes) == 0L) {
    return(invisible(NULL))
    }

    N <- stats::median(sizes)
    fraction <- n / N
    if (!is.finite(fraction) || fraction <= threshold) {
    return(invisible(NULL))
    }

    opening <- if (identical(when, "results")) {
    "These results come from a foreground of"
    } else {
    "Scanning a foreground of"
    }
    warning(
    opening, " ", n, " sequences against a background of ", N,
    " promoters (f = ", format(round(fraction, 3), nsmall = 2), "). ",
    "When the foreground size exceeds ", round(100 * threshold), "% of the ",
    "background (here ", max(1, round(threshold * N)), " sequences), the ",
    "z-score statistic Pscan relies on becomes increasingly less reliable. ",
    "See ?pscan, section \"Foreground size\", for why. ",
    "Set options(PscanR.foreground.max_fraction = ) to change this ",
    "threshold, or Inf to silence it.",
    call. = FALSE
    )
    invisible(NULL)
}

# `ps_fg_size` and `ps_bg_size` are declared "integer" in AllClasses.R, which
# constrains the type but not the length: `integer(0)` is a valid slot value,
# and a hand-built or previously serialised PSMatrix can carry a double. A
# vapply(..., integer(1L)) over them therefore turns a diagnostic into a hard
# error on objects that used to work, so anything that is not a single
# non-negative integer is reported here as NA_integer_ and left to the caller
# to ignore.
.ps_sizes <- function(pfms, accessor) {
    vapply(pfms, function(m) {
    size <- accessor(m)
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
        size < 0 || size > .Machine$integer.max) {
        return(NA_integer_)
    }
    as.integer(size)
    }, integer(1L))
}

.ps_check_short_matrix_file <- function(path) {
    if (file.access(path, mode = 4) != 0) {
    stop(sprintf("Cannot access file path: %s", path))
    }
    first_line <- readLines(path, n = 1)
    if (first_line != "[SHORT TFBS MATRIX]") {
    stop(sprintf("%s does not look like a Pscan .short_matrix file", path))
    }
}

.ps_check_bg_table <- function(x) {
    req_cols <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
    if (!all(req_cols %in% colnames(x))) {
    stop(
        "x does not contain required columns",
        "\"BG_SIZE\", \"BG_MEAN\", \"BG_STDEV\""
    )
    }
    if (!all(
    is.numeric(x$BG_SIZE), is.numeric(x$BG_MEAN),
    is.numeric(x$BG_STDEV)
    )) {
    stop("Required columns of x must be of numeric type")
    }
}

.ps_checks <- function(x, pfms, type) {
    # .ps_required_packages()

    if (type == 4 && !is(pfms, "PSMatrixList")) {
    stop("pfms is not an object of PSMatrixList class")
    }

    if (type == 4) {
    .ps_warn_missing_bg(pfms)
    }

    if (!is(pfms, "PFMatrixList") && !is(pfms, "PSMatrixList")) {
    stop("pfms is not an object of PFMatrixList or PSMatrixList class")
    }

    if (is.character(x) && type == 2) {
    .ps_check_short_matrix_file(x)
    } else if (!is(x, "DNAStringSet") && (type == 1 || type == 4)) {
    stop("x is not an object of DNAStringSet class")
    } else if (is.data.frame(x) && type == 3) {
    .ps_check_bg_table(x)
    }
}

#' @keywords internal
.ps_checks2 <- function(pfms, file = NULL, ...) {
    # .ps_required_packages()

    if (!is(pfms, "PSMatrixList")) {
    stop("pfms is not an object of PSMatrixList class")
    }

    if (!is.null(file) && !is.character(file)) {
    stop("file must be a character string indicating the file path")

    file_dir <- dirname(file)

    if (!file.exists(file)) {
        if (file.access(file_dir, mode = 2) != 0) {
        stop(sprintf("Can't write to directory: %s", file_dir))
        }
    } else {
        if (file.access(file, mode = 2) != 0) {
        stop(sprintf("Can't write to existing file: %s", file))
        }
    }
    }
}

#' @keywords internal
.clean_sequence <- function(x) {
    seq_widths <- Biostrings::width(x)
    ref_width <- max(Biostrings::width(x))

    diff_length_seq <- x[seq_widths != ref_width]

    if (length(diff_length_seq) != 0) {
    warning(paste(
        length(diff_length_seq), "sequences found with length
                    different from the reference. Removing the following
                    sequences:",
        paste(names(diff_length_seq), collapse = ", ")
    ))
    }
    x <- x[seq_widths == ref_width]

    freqs <- Biostrings::alphabetFrequency(x, as.prob = FALSE)
    n_proportions <- freqs[, "N"] / ref_width
    rem_names <- names(x[n_proportions > 0.5])

    if (length(rem_names) > 0) {
    warning(paste(
        "Found", length(rem_names), "sequences with more than 50% of N.
                    Removing the following sequences:",
        paste(rem_names, collapse = ", ")
    ))
    }
    x <- x[n_proportions <= 0.5]
    return(x)
}

# The legend is a slot of the list, not of its matrices, so every
# do.call(PSMatrixList, ...) rebuild loses it unless it is carried across by
# hand. Returns character() for a PFMatrixList, which has no such slot.
#' @keywords internal
.ps_legend_of <- function(x) {
    if (is(x, "PSMatrixList")) x@transcriptIDLegend else character()
}

# A full background stores the hits of every promoter in the universe. An
# ordinary background stores only their summary statistics, and a scan result
# stores foreground hits, so neither can serve a retrieval request.
#' @keywords internal
.ps_has_bg_scan <- function(pfms) {
    length(pfms) > 0 && all(vapply(
    pfms, function(m) length(ps_hits_score_bg(m)) > 0, logical(1)
    ))
}

# ps_bg_size() is set to the number of hits retained when a background is
# scanned (.ps_add_hits(), under `if (BG)`), and fullBG copies those same
# vectors into the _bg slots, so length(ps_hits_score_bg(x)) == ps_bg_size(x)
# by construction. A table that disagrees was computed on a different set of
# promoters, and applying it would leave the object with statistics from one
# universe and hits from another -- which pscan_fullBG() would then read as a
# background and answer from, silently.
#
# The slot is read directly rather than through ps_hits_score_bg(), which
# attaches the sequence names to every score on the way out: that is a copy per
# matrix for a length, and it would raise a names<- error of its own on an
# object already malformed enough to reach this check.
#' @keywords internal
.ps_check_bg_scan_size <- function(x, pfms) {
    ids <- vapply(pfms, ID, character(1L))
    stored <- vapply(pfms, function(m) length(m@ps_hits_score_bg), integer(1L))

    # Only matrices the table actually names, as .ps_bg_from_table() does, and
    # only those carrying a scan -- an ordinary PFMatrixList has none.
    compare <- ids %in% row.names(x) & stored > 0L
    if (!any(compare)) {
    return(invisible(NULL))
    }

    ids <- ids[compare]
    stored <- stored[compare]
    tabled <- as.integer(x[ids, "BG_SIZE"])

    bad <- which(!is.na(tabled) & tabled != stored)
    if (length(bad) == 0L) {
    return(invisible(NULL))
    }

    first <- bad[1L]
    stop(
    "Background table disagrees with the stored background scan for ",
    length(bad), " of ", length(ids), " matrices (e.g. ", ids[first],
    ": the table says ", tabled[first], " promoters, the stored scan holds ",
    stored[first], "). The two are equal for a background built by ",
    "ps_build_bg(), so this table was computed on a different set of ",
    "promoters from the one these matrices were scanned against.",
    call. = FALSE
    )
}

#' @keywords internal
.mapping_unique_names <- function(x, pfms) {
    original_names <- names(x)

    all_sequences_ID <- setNames(character(length(x)), original_names)

    # duplicate elimination
    unique_x <- BiocGenerics::unique(x)
    unique_names <- names(unique_x)

    x_char <- as.character(x)
    unique_x_char <- as.character(unique_x)

    all_sequences_ID <- setNames(
    unique_names[match(x_char, unique_x_char)],
    original_names
    )

    x <- .clean_sequence(x)
    removed_sequences <- setdiff(original_names, names(x))

    # NA corresponds to sequences eliminated by .clean_sequence()
    all_sequences_ID[removed_sequences] <- NA

    pfms@transcriptIDLegend <- all_sequences_ID

    return(pfms)
}

#' @keywords internal
#' @importFrom utils download.file
.ps_background_repository <- "Federico77z/PscanR_backgrounds"
.ps_background_sources <- c("experimenthub", "zenodo", "github")
.ps_experimenthub_package <- "PscanRBackgrounds"
# The archive backends distribute exactly one background version. The Hub
# title, the Zenodo file name, the checksum and the directory prefix inside
# the archive all describe that same release, so they are derived from a
# single constant and cross-checked in .ps_extract_background().
.ps_background_archive_version <- 2L
.ps_experimenthub_title <- paste0(
    "PscanR_backgrounds_v", .ps_background_archive_version
)
.ps_zenodo_record_id <- "21821764"
.ps_zenodo_archive_name <- paste0(.ps_experimenthub_title, ".zip")
.ps_zenodo_archive_sha256 <-
    "668b80839f2b81c7f48d11d0640a06fd58774c6cc9702ec73f4f54ebe9d0b625"

.ps_match_background_source <- function(source) {
    if (!is.character(source) || !length(source) || anyNA(source)) {
        stop(
            "source must be one of: ",
            paste(.ps_background_sources, collapse = ", "), call. = FALSE
        )
    }
    match.arg(tolower(source), .ps_background_sources)
}

.ps_background_raw_url <- function(path) {
    base <- getOption(
        "PscanR.background.base_url",
        paste0(
            "https://raw.githubusercontent.com/",
            .ps_background_repository,
            "/refs/heads/main"
        )
    )
    paste0(sub("/$", "", base), "/", sub("^/", "", path))
}

.ps_catalog_source <- function() {
    getOption(
        "PscanR.background.catalog",
        .ps_background_raw_url("catalog.tsv")
    )
}

.ps_hub_catalog_source <- function() {
    source <- getOption("PscanR.background.hub_catalog", NULL)
    if (is.null(source)) {
        source <- system.file(
            "extdata", "PscanR_background_catalog_v2.tsv",
            package = "PscanR"
        )
    }
    if (!length(source) || !nzchar(source)) {
        stop(
            "Bundled ExperimentHub background catalog is missing",
            call. = FALSE
        )
    }
    source
}

.ps_background_catalog <- function(source) {
    source <- .ps_match_background_source(source)
    catalog_source <- if (source == "github") {
        .ps_catalog_source()
    } else {
        .ps_hub_catalog_source()
    }
    .ps_read_bg_catalog(catalog_source)
}

.ps_read_bg_catalog <- function(source = .ps_catalog_source()) {
    required <- c(
        "status", "latest", "organism", "assembly", "upstream",
        "downstream", "jaspar_release", "background_version", "artifact",
        "artifact_sha256"
    )
    catalog <- tryCatch(
        utils::read.delim(
            source, stringsAsFactors = FALSE, check.names = FALSE,
            colClasses = "character", na.strings = character()
        ),
        error = function(error) {
            stop(
                "Unable to retrieve the PscanR background catalog: ",
                conditionMessage(error), call. = FALSE
            )
        }
    )
    missing <- setdiff(required, names(catalog))
    if (length(missing)) {
        stop(
            "Invalid PscanR background catalog; missing columns: ",
            paste(missing, collapse = ", "), call. = FALSE
        )
    }
    for (column in c(
        "upstream", "downstream", "jaspar_release", "background_version"
    )) {
        values <- catalog[[column]]
        if (any(!grepl("^[0-9]+$", values))) {
            stop(
                "Invalid integer values in catalog column '", column, "'",
                call. = FALSE
            )
        }
        catalog[[column]] <- as.integer(values)
    }
    catalog$latest <- tolower(catalog$latest) == "true"
    catalog
}

.ps_bg_target <- function(JASPAR_matrix, org, prom_reg, assembly, version) {
    if (!is.character(JASPAR_matrix) || length(JASPAR_matrix) != 1L ||
        is.na(JASPAR_matrix)) {
        stop("JASPAR_matrix must be a single string", call. = FALSE)
    }
    JASPAR_matrix <- toupper(JASPAR_matrix)
    if (!grepl("^JASPAR[0-9]{4}$", JASPAR_matrix)) {
        stop(
            "JASPAR_matrix must have the form JASPAR followed by a year",
            call. = FALSE
        )
    }
    if (!is.character(org) || length(org) != 1L || is.na(org) ||
        !org %in% c("hs", "mm", "at", "sc", "dm")) {
        stop("Unsupported organism code: ", org, call. = FALSE)
    }
    resolved_assembly <- switch(
        org, hs = assembly, mm = assembly, at = "TAIR9",
        sc = "sacCer3", dm = "dm6"
    )
    if (length(resolved_assembly) != 1L || is.na(resolved_assembly) ||
        !nzchar(resolved_assembly)) {
        stop("assembly is required for organism '", org, "'", call. = FALSE)
    }
    if (!is.numeric(prom_reg) || length(prom_reg) != 2L ||
        anyNA(prom_reg) || prom_reg[[1]] > 0 || prom_reg[[2]] < 0) {
        stop("prom_reg must be c(-upstream, downstream)", call. = FALSE)
    }
    version <- as.character(version)
    if (length(version) != 1L || is.na(version)) {
        stop("version must be 'latest' or a positive integer", call. = FALSE)
    }
    latest <- identical(tolower(version), "latest")
    if (!latest && !grepl("^[1-9][0-9]*$", version)) {
        stop("version must be 'latest' or a positive integer", call. = FALSE)
    }
    list(
        release = as.integer(sub("^JASPAR", "", JASPAR_matrix)),
        assembly = resolved_assembly,
        upstream = abs(as.integer(prom_reg[[1]])),
        downstream = abs(as.integer(prom_reg[[2]])),
        latest = latest,
        version = if (latest) NA_integer_ else as.integer(version),
        version_label = version
    )
}

.ps_resolve_bg_catalog <- function(
    catalog, JASPAR_matrix, org, prom_reg, assembly, version
) {
    target <- .ps_bg_target(
        JASPAR_matrix, org, prom_reg, assembly, version
    )
    candidates <- catalog[
        catalog$status == "validated" & catalog$organism == org &
            catalog$assembly == target$assembly &
            catalog$upstream == target$upstream &
            catalog$downstream == target$downstream &
            catalog$jaspar_release == target$release,
        , drop = FALSE
    ]
    if (target$latest) {
        candidates <- candidates[candidates$latest, , drop = FALSE]
    } else {
        candidates <- candidates[
            candidates$background_version == target$version, , drop = FALSE
        ]
    }
    if (nrow(candidates) != 1L) {
        stop(
            "No unique validated background matches JASPAR", target$release,
            ", ", target$assembly, ", ", target$upstream, "u_",
            target$downstream, "d, version ", target$version_label,
            ". Run get_availableBG(details = TRUE) to inspect the catalog.",
            call. = FALSE
        )
    }
    candidates[1L, , drop = FALSE]
}

.ps_verify_sha256 <- function(path, sha256, label = basename(path)) {
    if (is.null(sha256) || anyNA(sha256) || !nzchar(sha256)) {
        return(invisible(path))
    }
    actual <- unname(tools::sha256sum(path))
    if (!identical(tolower(actual), tolower(sha256))) {
        stop(label, " failed SHA-256 validation", call. = FALSE)
    }
    invisible(path)
}

.download_background <- function(file, destfile = NULL, sha256 = NULL) {
    if (is.null(destfile)) {
        destfile <- file.path(tempdir(), basename(file))
    }
    URL <- .ps_background_raw_url(file)
    # Download to a staging file so that a failed transfer or a checksum
    # mismatch never overwrites or removes an existing destfile.
    staged <- tempfile("PscanR-background-")
    on.exit(unlink(staged), add = TRUE)
    utils::download.file(URL, staged, mode = "wb", quiet = TRUE)
    .ps_verify_sha256(
        staged, sha256,
        label = paste("Downloaded background", basename(file))
    )
    if (!file.copy(staged, destfile, overwrite = TRUE)) {
        stop(
            "Could not save the downloaded background to ", destfile,
            call. = FALSE
        )
    }
    return(destfile)
}

.ps_zenodo_archive_url <- function() {
    getOption(
        "PscanR.background.zenodo_url",
        paste0(
            "https://zenodo.org/api/records/", .ps_zenodo_record_id,
            "/files/", .ps_zenodo_archive_name, "/content"
        )
    )
}

.ps_expected_archive_sha256 <- function() {
    getOption(
        "PscanR.background.archive_sha256", .ps_zenodo_archive_sha256
    )
}

.ps_background_cache <- function() {
    path <- getOption(
        "PscanR.background.cache",
        tools::R_user_dir("PscanR", which = "cache")
    )
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, mustWork = TRUE)
}

.ps_open_experimenthub <- function() {
    ExperimentHub::ExperimentHub()
}

.ps_fetch_experimenthub_archive <- function() {
    hub <- .ps_open_experimenthub()
    keep <- as.character(hub$preparerclass) == .ps_experimenthub_package &
        as.character(hub$title) == .ps_experimenthub_title
    if (sum(keep) != 1L) {
        stop(
            "ExperimentHub does not contain one unique ",
            .ps_experimenthub_title, " resource", call. = FALSE
        )
    }
    archive <- hub[[names(hub)[keep]]]
    if (!is.character(archive) || length(archive) != 1L ||
        !file.exists(archive)) {
        stop(
            "ExperimentHub returned an invalid background archive",
            call. = FALSE
        )
    }
    .ps_verify_sha256(
        archive, .ps_expected_archive_sha256(),
        label = "ExperimentHub background archive"
    )
    archive
}

.ps_fetch_zenodo_archive <- function() {
    cache <- BiocFileCache::BiocFileCache(
        cache = .ps_background_cache(), ask = FALSE
    )
    URL <- .ps_zenodo_archive_url()
    archive <- unname(BiocFileCache::bfcrpath(
        cache, rnames = URL, exact = TRUE
    ))
    if (length(archive) != 1L || !file.exists(archive)) {
        stop("Zenodo did not return a background archive", call. = FALSE)
    }
    tryCatch(
        .ps_verify_sha256(
            archive, .ps_expected_archive_sha256(),
            label = "Zenodo background archive"
        ),
        error = function(error) {
            record <- BiocFileCache::bfcquery(
                cache, URL, field = "rname", exact = TRUE
            )
            if (nrow(record)) BiocFileCache::bfcremove(cache, record$rid)
            stop(conditionMessage(error), call. = FALSE)
        }
    )
    archive
}

.ps_fetch_background_archive <- function(source) {
    source <- .ps_match_background_source(source)
    if (source == "zenodo") return(.ps_fetch_zenodo_archive())
    if (source == "github") {
        stop(
            "GitHub backgrounds are not stored in the Zenodo archive",
            call. = FALSE
        )
    }
    tryCatch(
        .ps_fetch_experimenthub_archive(),
        error = function(error) {
            warning(
                "ExperimentHub background retrieval failed: ",
                conditionMessage(error),
                ". Using the immutable Zenodo fallback.", call. = FALSE
            )
            .ps_fetch_zenodo_archive()
        }
    )
}

.ps_extract_background <- function(archive, entry, destfile = NULL) {
    version <- as.integer(entry$background_version[[1]])
    if (!identical(version, .ps_background_archive_version)) {
        stop(
            "Background version ", version, " is not distributed in the ",
            .ps_zenodo_archive_name, " archive, which holds version ",
            .ps_background_archive_version,
            ". Use source = \"github\" for other versions.", call. = FALSE
        )
    }
    member <- paste0(
        .ps_experimenthub_title, "/", sub("^/", "", entry$artifact[[1]])
    )
    members <- utils::unzip(archive, list = TRUE)$Name
    if (sum(members == member) != 1L) {
        stop(
            "Background archive does not contain one unique ", member,
            call. = FALSE
        )
    }
    staging <- tempfile("PscanR-background-")
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    utils::unzip(archive, files = member, exdir = staging)
    extracted <- do.call(
        file.path,
        as.list(c(staging, strsplit(member, "/", fixed = TRUE)[[1]]))
    )
    # Validate in the staging directory: destfile may be an existing user file,
    # which must not be overwritten or removed on account of a bad member.
    .ps_verify_sha256(
        extracted, entry$artifact_sha256[[1]],
        label = paste("Background", basename(entry$artifact[[1]]))
    )
    if (is.null(destfile)) {
        destfile <- tempfile(fileext = ".txt")
    }
    if (!file.copy(extracted, destfile, overwrite = TRUE)) {
        stop(
            "Could not save the selected background to ", destfile,
            call. = FALSE
        )
    }
    destfile
}

.ps_retrieve_background <- function(entry, source, destfile = NULL) {
    source <- .ps_match_background_source(source)
    if (source == "github") {
        return(.download_background(
            file = entry$artifact[[1]], destfile = destfile,
            sha256 = entry$artifact_sha256[[1]]
        ))
    }
    archive <- .ps_fetch_background_archive(source)
    .ps_extract_background(archive, entry, destfile)
}

.check_seq_duplicated <- function(x) {
    for (i in seq_along(x)) {
    if (names(x[i]) != x[i]) {
        warning(sprintf(
        "%s will be evaluated instead of %s since they have the same %s",
        x[i], names(x[i]), "promoter region"
        ))
    }
    }
}

# .ps_required_packages <- function()
# {
#  require("Biostrings")
#  require("TFBSTools")
#  require("BiocParallel")
#  require("BSDA")
# }
