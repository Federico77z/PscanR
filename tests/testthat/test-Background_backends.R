make_background_archive_fixture <- function(background, version = 2L) {
    root <- tempfile("pscan-background-archive-")
    release <- file.path(root, paste0("PscanR_backgrounds_v", version))
    dir.create(file.path(release, "backgrounds"), recursive = TRUE)
    target <- file.path(release, "backgrounds", basename(background))
    stopifnot(file.copy(background, target))
    archive <- tempfile(fileext = ".zip")
    old <- setwd(root)
    on.exit(setwd(old), add = TRUE)
    status <- utils::zip(
        archive,
        files = file.path(
            basename(release), "backgrounds", basename(background)
        ),
        flags = "-X -q"
    )
    setwd(old)
    stopifnot(identical(status, 0L))
    archive
}

test_that("bundled Hub catalog contains the complete version-2 release", {
    details <- get_availableBG(details = TRUE)
    expect_identical(nrow(details), 105L)
    expect_true(all(details$status == "validated"))
    expect_true(all(details$latest))
    expect_true(all(details$background_version == 2L))
    expect_error(
        PscanR:::.ps_match_background_source("automatic"),
        "should be one of"
    )
})

test_that("ExperimentHub failure falls back to Zenodo with a warning", {
    archive <- tempfile(fileext = ".zip")
    writeBin(charToRaw("fixture"), archive)
    local_mocked_bindings(
        .ps_fetch_experimenthub_archive = function() stop("Hub unavailable"),
        .ps_fetch_zenodo_archive = function() archive,
        .package = "PscanR"
    )

    expect_warning(
        result <- PscanR:::.ps_fetch_background_archive("experimenthub"),
        "Using the immutable Zenodo fallback"
    )
    expect_identical(result, archive)
})

test_that("archive and GitHub backends yield identical background objects", {
    background <- system.file(
        "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
        package = "PscanR"
    )
    matrices <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
    archive <- make_background_archive_fixture(background)
    entry <- data.frame(
        background_version = 2L,
        artifact = file.path("backgrounds", basename(background)),
        artifact_sha256 = unname(tools::sha256sum(background)),
        stringsAsFactors = FALSE
    )

    archive_path <- PscanR:::.ps_extract_background(archive, entry)

    github_root <- tempfile("pscan-github-fixture-")
    dir.create(file.path(github_root, "BG_files"), recursive = TRUE)
    github_source <- file.path(github_root, "BG_files", basename(background))
    expect_true(file.copy(background, github_source))
    old <- options(
        PscanR.background.base_url = paste0(
            "file://", normalizePath(github_root)
        )
    )
    on.exit(options(old), add = TRUE)
    github_path <- PscanR:::.download_background(
        file.path("BG_files", basename(background)),
        sha256 = entry$artifact_sha256
    )

    archive_result <- ps_retrieve_bg_from_file(archive_path, matrices)
    github_result <- ps_retrieve_bg_from_file(github_path, matrices)
    expect_identical(archive_result, github_result)
})

test_that("archive extraction validates member checksums and destfile", {
    background <- system.file(
        "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
        package = "PscanR"
    )
    archive <- make_background_archive_fixture(background)
    entry <- data.frame(
        background_version = 2L,
        artifact = file.path("backgrounds", basename(background)),
        artifact_sha256 = unname(tools::sha256sum(background)),
        stringsAsFactors = FALSE
    )
    destination <- tempfile(fileext = ".txt")
    result <- PscanR:::.ps_extract_background(archive, entry, destination)
    expect_identical(result, destination)
    expect_identical(readLines(result), readLines(background))

    entry$artifact_sha256 <- paste(rep("0", 64L), collapse = "")
    bad_destination <- tempfile(fileext = ".txt")
    expect_error(
        PscanR:::.ps_extract_background(archive, entry, bad_destination),
        "SHA-256"
    )
    expect_false(file.exists(bad_destination))

    # A destfile that already exists must survive a checksum failure intact.
    preserved <- tempfile(fileext = ".txt")
    writeLines("existing user content", preserved)
    expect_error(
        PscanR:::.ps_extract_background(archive, entry, preserved), "SHA-256"
    )
    expect_true(file.exists(preserved))
    expect_identical(readLines(preserved), "existing user content")
})

test_that("archive extraction rejects versions the archive does not carry", {
    background <- system.file(
        "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
        package = "PscanR"
    )
    archive <- make_background_archive_fixture(background)
    entry <- data.frame(
        background_version = 1L,
        artifact = file.path("backgrounds", basename(background)),
        artifact_sha256 = unname(tools::sha256sum(background)),
        stringsAsFactors = FALSE
    )
    expect_error(
        PscanR:::.ps_extract_background(archive, entry),
        "is not distributed in the"
    )
})

test_that("version 1 is available only through the explicit GitHub catalog", {
    hub_catalog <- get_availableBG(details = TRUE, source = "experimenthub")
    expect_false(any(hub_catalog$background_version == 1L))

    catalog_path <- tempfile(fileext = ".tsv")
    catalog <- data.frame(
        status = c("validated", "validated"), latest = c(FALSE, TRUE),
        organism = c("hs", "hs"), assembly = c("hg38", "hg38"),
        upstream = c(950L, 950L), downstream = c(50L, 50L),
        jaspar_release = c(2024L, 2024L),
        background_version = c(1L, 2L),
        artifact = c(
            "BG_files/J2024_hg38_950u_50d_UCSC.psbg1.txt",
            "BG_files/J2024_hg38_950u_50d_UCSC.psbg2.txt"
        ),
        artifact_sha256 = c("one", "two"), stringsAsFactors = FALSE
    )
    utils::write.table(
        catalog, catalog_path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    old <- options(PscanR.background.catalog = catalog_path)
    on.exit(options(old), add = TRUE)
    github_catalog <- get_availableBG(details = TRUE, source = "github")
    expect_true(any(github_catalog$background_version == 1L))
})
