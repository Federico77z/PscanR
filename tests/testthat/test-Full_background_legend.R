## The transcriptIDLegend is what makes a PSMatrixList a full background. It
## lives on the list, not on its matrices, so every rebuild has to carry it
## across by hand. These tests pin that it does.

full_pfms <- function() {
    p <- system.file("extdata", "full_pfms.rds", package = "PscanR")
    skip_if_not(nzchar(p), "full_pfms.rds example data not installed")
    readRDS(p)
}

prom_seq <- function() {
    p <- system.file("extdata", "prom_seq.rds", package = "PscanR")
    skip_if_not(nzchar(p), "prom_seq.rds example data not installed")
    readRDS(p)
}

# The bundled background is a 36-promoter toy, so any usable foreground is a
# large share of it. The fraction warning is not what these tests are about.
silence_fraction <- function() options(PscanR.foreground.max_fraction = Inf)

test_that("PSMatrixList() honours its transcriptIDLegend argument", {
    full <- full_pfms()
    legend <- c(ID1 = "ID2", ID2 = "ID2", ID3 = "ID2")

    out <- PSMatrixList(full[[1]], full[[2]], transcriptIDLegend = legend)

    expect_identical(transcriptIDLegend(out), legend)
    # The default still marks an ordinary list.
    expect_identical(transcriptIDLegend(PSMatrixList(full[[1]])), character())
    expect_error(
        PSMatrixList(full[[1]], transcriptIDLegend = 1:3),
        "must be a character vector"
    )
})

test_that("a full background survives the text round trip with its legend", {
    full <- full_pfms()
    old <- silence_fraction()
    on.exit(options(old), add = TRUE)

    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f), add = TRUE)
    ps_write_bg_to_file(full, f)

    back <- ps_retrieve_bg_from_file(f, full)

    expect_identical(transcriptIDLegend(back), transcriptIDLegend(full))
    # The file carries only the three summary statistics; the per-promoter
    # scan is still the one the matrices came in with.
    expect_identical(
        ps_hits_score_bg(back[[1]]), ps_hits_score_bg(full[[1]])
    )

    ids <- utils::head(names(transcriptIDLegend(full)), 10)
    from_back <- suppressWarnings(pscan_fullBG(ids, back, quiet = TRUE))
    from_full <- suppressWarnings(pscan_fullBG(ids, full, quiet = TRUE))
    expect_identical(
        ps_hits_score(from_back[[1]]), ps_hits_score(from_full[[1]])
    )
})

test_that("pscan, pscan_fullBG and PscanFiltered preserve the legend", {
    full <- full_pfms()
    prom <- prom_seq()
    legend <- transcriptIDLegend(full)
    old <- silence_fraction()
    on.exit(options(old), add = TRUE)

    scanned <- pscan(prom[seq_len(5)], full[seq_len(2)],
                     BPPARAM = BiocParallel::SerialParam())
    expect_identical(transcriptIDLegend(scanned), legend)

    ids <- utils::head(names(legend), 10)
    retrieved <- suppressWarnings(pscan_fullBG(ids, full, quiet = TRUE))
    expect_identical(transcriptIDLegend(retrieved), legend)

    filtered <- PscanFiltered(prom[seq_len(20)], full[[1]], n = 1,
                              background = full[seq_len(2)],
                              BPPARAM = BiocParallel::SerialParam())
    expect_identical(transcriptIDLegend(filtered), legend)
})

test_that("pscan_fullBG rejects a legend without a background scan", {
    full <- full_pfms()
    prom <- prom_seq()
    old <- silence_fraction()
    on.exit(options(old), add = TRUE)

    # A pscan() result now carries the background's legend, so the legend on
    # its own no longer identifies a full background. Reading one as a
    # background would answer with foreground hits over a subset.
    scanned <- pscan(prom[seq_len(5)], full[seq_len(2)],
                     BPPARAM = BiocParallel::SerialParam())
    expect_true(length(transcriptIDLegend(scanned)) > 0)
    expect_error(
        pscan_fullBG(names(transcriptIDLegend(full))[1], scanned, quiet = TRUE),
        "no per-promoter background scan"
    )

    # The legend is still required, and says so separately.
    legendless <- full
    legendless@transcriptIDLegend <- character()
    expect_error(
        pscan_fullBG(names(transcriptIDLegend(full))[1], legendless,
                     quiet = TRUE),
        "empty transcriptIDLegend"
    )

    expect_error(pscan_fullBG("NM_000546", full[[1]]), "must be a PSMatrixList")
})

test_that("ps_build_bg still sets the legend only for a full background", {
    prom <- prom_seq()
    full <- full_pfms()
    x <- prom[seq_len(8)]

    ordinary <- ps_build_bg(x, full[seq_len(2)],
                            BPPARAM = BiocParallel::SerialParam())
    expect_identical(transcriptIDLegend(ordinary), character())

    complete <- ps_build_bg(x, full[seq_len(2)],
                            BPPARAM = BiocParallel::SerialParam(),
                            fullBG = TRUE)
    legend <- transcriptIDLegend(complete)
    expect_identical(names(legend), names(x))
    # Each name maps to the representative kept by unique(); here every
    # promoter is distinct, so each maps to itself.
    expect_identical(unname(legend), names(x))
})
