## Tests for background construction, import/export and the ExperimentHub-backed
## retrieval helpers.

j2020 <- function() {
    p <- system.file("extdata", "J2020.rds", package = "PscanR")
    skip_if_not(nzchar(p), "J2020.rds example data not installed")
    readRDS(p)
}

test_that("ps_retrieve_bg_from_file loads bundled background statistics", {
    bg_file <- system.file(
        "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", package = "PscanR")
    skip_if_not(nzchar(bg_file), "bundled background file not installed")
    J2020 <- j2020()

    bg <- ps_retrieve_bg_from_file(bg_file, J2020)

    expect_s4_class(bg, "PSMatrixList")
    expect_equal(length(bg), length(J2020))

    avg <- ps_bg_avg(bg[[1]])
    std <- ps_bg_std_dev(bg[[1]])
    size <- ps_bg_size(bg[[1]])
    expect_true(is.finite(avg) && avg > 0 && avg < 1)
    expect_true(is.finite(std) && std > 0)
    expect_true(size > 0L)

    # A non-existent file is rejected
    expect_error(ps_retrieve_bg_from_file(tempfile(), J2020[1:2]))
})

test_that("ps_build_bg_from_table assigns statistics and round-trips", {
    J2020 <- j2020()
    J2020_subset <- J2020[1:3]
    ids <- c("MA0004.1", "MA0006.1", "MA0019.1")

    background_data <- data.frame(
        BG_SIZE = c(500L, 450L, 480L),
        BG_MEAN = c(0.30, 0.25, 0.35),
        BG_STDEV = c(0.05, 0.07, 0.06),
        row.names = ids
    )

    bg <- ps_build_bg_from_table(background_data, J2020_subset)
    expect_s4_class(bg, "PSMatrixList")
    expect_equal(length(bg), 3L)
    expect_equal(ps_bg_size(bg[[1]]), 500L)
    expect_equal(ps_bg_avg(bg[[2]]), 0.25)
    expect_equal(ps_bg_std_dev(bg[[3]]), 0.06)

    # ps_get_bg_table is the inverse of ps_build_bg_from_table
    tab <- ps_get_bg_table(bg)
    expect_s3_class(tab, "data.frame")
    expect_equal(colnames(tab), c("BG_SIZE", "BG_MEAN", "BG_STDEV"))
    expect_equal(tab[ids, "BG_MEAN"], background_data[ids, "BG_MEAN"])

    # Row/PFM count mismatch only warns (does not error)
    expect_warning(ps_build_bg_from_table(background_data, J2020_subset[1:2]))
})

test_that("ps_write_bg_to_file and ps_retrieve_bg_from_file round-trip", {
    J2020 <- j2020()
    sub <- J2020[1:2]

    PSM1 <- PSMatrix(pfm = sub[[1]], ps_bg_avg = 0.25, ps_fg_avg = 0.5,
                     ps_bg_std_dev = 0.05, ps_bg_size = 250L)
    PSM2 <- PSMatrix(pfm = sub[[2]], ps_bg_avg = 0.30, ps_fg_avg = 0.5,
                     ps_bg_std_dev = 0.06, ps_bg_size = 300L)
    psl <- PSMatrixList(PSM1, PSM2)
    names(psl) <- c(ID(PSM1), ID(PSM2))

    f <- tempfile(fileext = ".txt")
    on.exit(unlink(f), add = TRUE)
    expect_invisible(ps_write_bg_to_file(psl, f))

    # The file carries the expected header
    expect_identical(readLines(f, n = 1), "[SHORT TFBS MATRIX]")

    back <- ps_retrieve_bg_from_file(f, sub)
    expect_s4_class(back, "PSMatrixList")
    expect_equal(ps_bg_avg(back[[1]]), 0.25)
    expect_equal(ps_bg_size(back[[2]]), 300L)
})

test_that("get_availableBG queries ExperimentHub", {
    skip_if_offline()
    skip_if_not_installed("ExperimentHub")
    bg <- tryCatch(get_availableBG(), error = function(e) skip(conditionMessage(e)))
    expect_type(bg, "character")
    # Resources, once published, follow the J<year>_<assembly>_..._.psbg<v>.txt
    # naming convention; before publication the result may legitimately be empty.
    if (length(bg) > 0L) {
        expect_true(all(grepl("\\.txt$", bg)))
    }
})
