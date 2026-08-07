make_catalog_fixture <- function(path) {
    catalog <- data.frame(
        schema_version = c(1L, 1L),
        status = c("validated", "validated"),
        latest = c(FALSE, TRUE),
        organism = c("hs", "hs"),
        assembly = c("hg38", "hg38"),
        upstream = c(950L, 950L),
        downstream = c(50L, 50L),
        jaspar_release = c(2024L, 2024L),
        tax_group = c("vertebrates", "vertebrates"),
        background_version = c(1L, 2L),
        artifact = c(
            "BG_files/J2024_hg38_950u_50d_UCSC.psbg1.txt",
            "BG_files/J2024_hg38_950u_50d_UCSC.psbg2.txt"
        ),
        artifact_sha256 = c("one", "two"),
        stringsAsFactors = FALSE
    )
    utils::write.table(
        catalog, path, sep = "\t", quote = FALSE, row.names = FALSE
    )
    catalog
}

test_that("background catalog supports legacy and detailed listings", {
    catalog_path <- tempfile(fileext = ".tsv")
    make_catalog_fixture(catalog_path)
    old <- options(PscanR.background.catalog = catalog_path)
    on.exit(options(old), add = TRUE)

    files <- get_availableBG(source = "github")
    expect_identical(length(files), 2L)
    expect_identical(
        files[[2]], "J2024_hg38_950u_50d_UCSC.psbg2.txt"
    )
    details <- get_availableBG("psbg2", details = TRUE, source = "github")
    expect_s3_class(details, "data.frame")
    expect_identical(nrow(details), 1L)
    expect_true(details$latest)
    expect_error(
        get_availableBG("missing", source = "github"), "Found 0 matches"
    )
    expect_error(
        get_availableBG(details = NA, source = "github"), "details must be"
    )
})

test_that("background catalog resolves latest and pinned versions", {
    catalog_path <- tempfile(fileext = ".tsv")
    make_catalog_fixture(catalog_path)
    catalog <- PscanR:::.ps_read_bg_catalog(catalog_path)

    latest <- PscanR:::.ps_resolve_bg_catalog(
        catalog, "jaspar2024", "hs", c(-950, 50), "hg38", "latest"
    )
    pinned <- PscanR:::.ps_resolve_bg_catalog(
        catalog, "JASPAR2024", "hs", c(-950, 50), "hg38", "1"
    )
    expect_identical(latest$background_version, 2L)
    expect_identical(pinned$background_version, 1L)
    expect_error(
        PscanR:::.ps_resolve_bg_catalog(
            catalog, "JASPAR2018", "hs", c(-950, 50), "hg38", "latest"
        ),
        "No unique validated background"
    )
    expect_error(
        PscanR:::.ps_resolve_bg_catalog(
            catalog, "JASPAR2024", "hs", c(950, 50), "hg38", "latest"
        ),
        "prom_reg"
    )
    expect_error(
        PscanR:::.ps_resolve_bg_catalog(
            catalog, "JASPAR2024", "hs", c(-950, 50), "hg38", "zero"
        ),
        "version must be"
    )
})

test_that("background downloads enforce catalog checksums", {
    source_root <- tempfile("pscan-background-")
    dir.create(file.path(source_root, "BG_files"), recursive = TRUE)
    source_file <- file.path(source_root, "BG_files", "fixture.txt")
    writeLines(c("[SHORT TFBS MATRIX]", "MA0001.1\t10\t0.5\t0.1"), source_file)
    checksum <- unname(tools::sha256sum(source_file))
    base_url <- paste0("file://", normalizePath(source_root))
    old <- options(PscanR.background.base_url = base_url)
    on.exit(options(old), add = TRUE)

    destination <- tempfile(fileext = ".txt")
    result <- PscanR:::.download_background(
        "BG_files/fixture.txt", destination, checksum
    )
    expect_identical(result, destination)
    expect_identical(readLines(destination), readLines(source_file))

    corrupt_destination <- tempfile(fileext = ".txt")
    bad_checksum <- paste(rep("0", 64), collapse = "")
    expect_error(
        PscanR:::.download_background(
            "BG_files/fixture.txt", corrupt_destination, bad_checksum
        ),
        "SHA-256"
    )
    expect_false(file.exists(corrupt_destination))
})

test_that("malformed catalogs fail before network or matrix loading", {
    malformed <- tempfile(fileext = ".tsv")
    writeLines("status\tlatest", malformed)
    expect_error(
        PscanR:::.ps_read_bg_catalog(malformed),
        "missing columns"
    )
})
