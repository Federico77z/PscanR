test_that("argument checks stop with one wording and no internal call", {
    err <- function(expr) tryCatch(expr, error = identity)

    e <- err(.ps_assert_class(1, "DNAStringSet", "x"))
    expect_identical(conditionMessage(e), "'x' must be a DNAStringSet")
    expect_null(conditionCall(e))
    expect_error(.ps_assert_class(1, c("PSMatrix", "PSMatrixList"), "x"),
        "'x' must be a PSMatrix or PSMatrixList")
    seqs <- Biostrings::DNAStringSet("ACGT")
    expect_identical(.ps_assert_class(seqs, "DNAStringSet", "x"), seqs)

    expect_error(.ps_assert_character(1, "ID"),
        "'ID' must be a character vector")
    for (bad in list(1, c("a", "b"), NA_character_)) {
        expect_error(.ps_assert_string(bad, "keyword"),
            "'keyword' must be a single string")
    }
    for (bad in list("TRUE", c(TRUE, FALSE), NA)) {
        expect_error(.ps_assert_flag(bad, "details"),
            "'details' must be TRUE or FALSE")
    }
    expect_error(.ps_assert_number(2, "FDR", 0, 1),
        "'FDR' must be a single number between 0 and 1")
    expect_error(.ps_assert_number(0, "n", min = 1),
        "'n' must be a single number of at least 1")
    expect_error(.ps_assert_number("1", "n"), "'n' must be a single number$")
    expect_identical(.ps_assert_number(0.05, "FDR", 0, 1), 0.05)
})

test_that("an output file must be writable", {
    out <- tempfile(fileext = ".txt")
    expect_identical(.ps_assert_writable_file(out, "file"), out)
    missing_dir <- file.path(tempfile(), "sub", "bg.txt")
    expect_error(.ps_assert_writable_file(missing_dir, "file"),
        "Cannot write to the directory")
    expect_error(.ps_assert_writable_file(c("a", "b"), "file"),
        "'file' must be a single string")
})

test_that("exported functions check the arguments they used to let through", {
    J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
    expect_error(ps_retrieve_bg_from_file(42, J2020),
        "'file' must be a single string")
    expect_error(ps_build_bg_from_table(matrix(1), J2020),
        "'x' must be a data.frame")
    bg_file <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
        package = "PscanR")
    bg <- ps_retrieve_bg_from_file(bg_file, J2020)
    expect_error(
        ps_write_bg_to_file(bg, file.path(tempfile(), "sub", "bg.txt")),
        "Cannot write to the directory"
    )
})
