bundled_background <- function() {
    system.file(
        "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR"
    )
}

background_entry <- function(path = bundled_background(),
                             file = "backgrounds/J2020_hg38_200u_50d_UCSC.psbg2.txt") {
    data.frame(
        background_version = 2L,
        artifact = file,
        artifact_sha256 = unname(tools::sha256sum(path)),
        stringsAsFactors = FALSE
    )
}

mock_hub <- function(paths, titles, packages = rep(
    "PscanRBackgrounds", length(paths)
)) {
    structure(
        setNames(as.list(paths), paste0("EH", seq_along(paths))),
        class = "pscan_test_hub", title = titles, preparerclass = packages
    )
}

test_that("bundled catalog contains the complete version-2 release", {
    details <- ps_available_bg(details = TRUE)
    expect_identical(nrow(details), 105L)
    expect_true(all(details$status == "validated"))
    expect_true(all(details$latest))
    expect_true(all(details$background_version == 2L))
    expect_error(
        PscanR:::.ps_match_background_source("github"),
        "should be one of"
    )
})

make_background_archive <- function(background, file) {
    root <- tempfile("pscan-background-archive-")
    release <- file.path(root, "PscanR_backgrounds_v2", "backgrounds")
    dir.create(release, recursive = TRUE)
    stopifnot(file.copy(background, file.path(release, file)))
    archive <- tempfile(fileext = ".zip")
    old <- setwd(root)
    on.exit(setwd(old), add = TRUE)
    status <- utils::zip(
        archive, files = file.path("PscanR_backgrounds_v2", "backgrounds", file),
        flags = "-X -q"
    )
    stopifnot(identical(status, 0L))
    archive
}

test_that("catalog entries map to unique Hub titles", {
    expect_identical(
        PscanR:::.ps_hub_title(background_entry()),
        "PscanR_bg_v2_J2020_hg38_200u_50d_UCSC"
    )
    catalog <- ps_available_bg(details = TRUE)
    titles <- vapply(
        seq_len(nrow(catalog)),
        function(i) PscanR:::.ps_hub_title(catalog[i, ]),
        character(1L)
    )
    expect_false(anyDuplicated(titles) > 0L)
})

test_that("archive extraction validates the member and its version", {
    background <- bundled_background()
    entry <- background_entry()
    archive <- make_background_archive(background, basename(entry$artifact))
    extracted <- PscanR:::.ps_extract_background(archive, entry)
    expect_identical(readLines(extracted), readLines(background))

    bad <- entry
    bad$artifact_sha256 <- paste(rep("0", 64L), collapse = "")
    expect_error(PscanR:::.ps_extract_background(archive, bad), "SHA-256")
    missing <- entry
    missing$artifact <- "backgrounds/absent.psbg2.txt"
    expect_error(
        PscanR:::.ps_extract_background(archive, missing), "one unique"
    )
    old <- entry
    old$background_version <- 1L
    expect_error(
        PscanR:::.ps_extract_background(archive, old), "is not distributed"
    )
})

test_that("Hub lookup requires one matching resource and a valid file", {
    local_mocked_s3_method("$", "pscan_test_hub", function(x, name) {
        attr(x, name)
    })
    entry <- background_entry()
    title <- PscanR:::.ps_hub_title(entry)
    path <- bundled_background()
    hub <- mock_hub(path, title)
    local_mocked_bindings(
        .ps_open_experimenthub = function() hub, .package = "PscanR"
    )
    expect_identical(PscanR:::.ps_fetch_experimenthub_background(entry), path)
    hub <- mock_hub(path, "another-resource")
    expect_error(
        PscanR:::.ps_fetch_experimenthub_background(entry), "one unique"
    )
    hub <- mock_hub(rep(path, 2L), rep(title, 2L))
    expect_error(
        PscanR:::.ps_fetch_experimenthub_background(entry), "one unique"
    )
    hub <- mock_hub(rep(path, 2L), c(NA, title))
    expect_identical(PscanR:::.ps_fetch_experimenthub_background(entry), path)
    hub <- mock_hub(tempfile(), title)
    expect_error(
        PscanR:::.ps_fetch_experimenthub_background(entry), "invalid file"
    )
    hub <- mock_hub(NA_character_, title)
    expect_error(
        PscanR:::.ps_fetch_experimenthub_background(entry), "invalid file"
    )
    entry$artifact_sha256 <- paste(rep("0", 64L), collapse = "")
    hub <- mock_hub(path, title)
    expect_error(
        PscanR:::.ps_fetch_experimenthub_background(entry), "SHA-256"
    )
})

test_that("Zenodo archive is verified and a corrupt cache is evicted", {
    cache <- BiocFileCache::BiocFileCache(tempfile("pscan-cache-"), ask = FALSE)
    url <- "https://example.invalid/pscan-test.zip"
    corrupt <- tempfile(fileext = ".zip")
    writeBin(charToRaw("corrupt archive"), corrupt)
    BiocFileCache::bfcadd(cache, rname = url, fpath = corrupt, rtype = "local")
    local_mocked_bindings(
        .ps_background_cache = function() BiocFileCache::bfccache(cache),
        .ps_zenodo_archive_url = function() url,
        .package = "PscanR"
    )
    expect_error(PscanR:::.ps_fetch_zenodo_archive(), "SHA-256")
    expect_equal(nrow(BiocFileCache::bfcquery(cache, url, exact = TRUE)), 0L)
})

test_that("ExperimentHub failure falls back to Zenodo with a warning", {
    path <- bundled_background()
    local_mocked_bindings(
        .ps_fetch_experimenthub_background = function(entry) {
            stop("Hub unavailable")
        },
        .ps_fetch_zenodo_background = function(entry) path,
        .package = "PscanR"
    )
    expect_warning(
        result <- PscanR:::.ps_fetch_background(
            background_entry(), "experimenthub"
        ),
        "Using the Zenodo fallback"
    )
    expect_identical(result, path)
})

test_that("failure of both backends is reported", {
    local_mocked_bindings(
        .ps_fetch_experimenthub_background = function(entry) {
            stop("Hub unavailable")
        },
        .ps_fetch_zenodo_background = function(entry) {
            stop("Zenodo unavailable")
        },
        .package = "PscanR"
    )
    expect_warning(
        expect_error(
            PscanR:::.ps_fetch_background(background_entry(), "experimenthub"),
            "Zenodo unavailable"
        ),
        "Hub unavailable.*Zenodo fallback"
    )
})

test_that("a retrieved background is saved to destfile and read back", {
    path <- bundled_background()
    matrices <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
    local_mocked_bindings(
        .ps_fetch_zenodo_background = function(entry) path,
        .package = "PscanR"
    )
    destination <- tempfile(fileext = ".txt")
    result <- PscanR:::.ps_retrieve_background(
        background_entry(), "zenodo", destination
    )
    expect_identical(result, destination)
    expect_identical(readLines(destination), readLines(path))
    expect_identical(
        ps_retrieve_bg_from_file(destination, matrices),
        ps_retrieve_bg_from_file(path, matrices)
    )
    expect_identical(
        PscanR:::.ps_retrieve_background(background_entry(), "zenodo"), path
    )
})

test_that("an existing destfile survives a failed retrieval", {
    local_mocked_bindings(
        .ps_fetch_zenodo_background = function(entry) stop("SHA-256 mismatch"),
        .package = "PscanR"
    )
    preserved <- tempfile(fileext = ".txt")
    writeLines("existing user content", preserved)
    expect_error(
        PscanR:::.ps_retrieve_background(
            background_entry(), "zenodo", preserved
        ),
        "SHA-256"
    )
    expect_identical(readLines(preserved), "existing user content")
})
