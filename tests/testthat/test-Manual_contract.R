test_that("PSMatrix rejects unused constructor arguments", {
  motifs <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))

  expect_error(
    PSMatrix(motifs[[1]], misspelled_argument = 1),
    "unused argument.*misspelled_argument"
  )
})

test_that("PSMatrix validity covers scalar and hit-vector invariants", {
  foreground <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  full <- readRDS(system.file("extdata", "full_pfm1.rds", package = "PscanR"))

  zero_sd <- foreground
  zero_sd@ps_bg_std_dev <- 0
  expect_match(validPSMatrix(zero_sd), "Background stddev")

  bad_foreground <- foreground
  bad_foreground@ps_hits_oligo <- bad_foreground@ps_hits_oligo[-1]
  expect_match(validPSMatrix(bad_foreground), "foreground hit vectors")

  bad_background <- full
  bad_background@ps_hits_pos_bg <- bad_background@ps_hits_pos_bg[-1]
  expect_match(validPSMatrix(bad_background), "full-background hit vectors")

  bad_size <- full
  bad_size@ps_bg_size <- bad_size@ps_bg_size + 1L
  expect_match(validPSMatrix(bad_size), "does not equal ps_bg_size")
})

test_that("withDimnames controls accessor names", {
  foreground <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  full <- readRDS(system.file("extdata", "full_pfm1.rds", package = "PscanR"))

  foreground_accessors <- list(
    ps_hits_score, ps_hits_z, ps_hits_pos, ps_hits_strand, ps_hits_oligo
  )
  for (accessor in foreground_accessors) {
    expect_false(is.null(names(accessor(foreground))))
    expect_null(names(accessor(foreground, withDimnames = FALSE)))
  }

  background_accessors <- list(
    ps_hits_score_bg, ps_hits_pos_bg, ps_hits_strand_bg, ps_hits_oligo_bg
  )
  for (accessor in background_accessors) {
    expect_false(is.null(names(accessor(full))))
    expect_null(names(accessor(full, withDimnames = FALSE)))
  }

  named_table <- ps_hits_table(foreground)
  unnamed_table <- ps_hits_table(foreground, withDimnames = FALSE)
  expect_setequal(row.names(named_table), ps_seq_names(foreground))
  expect_identical(row.names(unnamed_table), as.character(seq_len(nrow(unnamed_table))))
  expect_type(named_table$OLIGO, "character")
})

test_that("reported coordinates assign pos_shift to the first base", {
  foreground <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))

  expect_identical(
    unname(ps_hits_pos(foreground)),
    foreground@ps_hits_pos - 1L
  )
  expect_identical(
    unname(ps_hits_pos(foreground, pos_shift = -200L)),
    foreground@ps_hits_pos - 1L - 200L
  )

  first_base <- foreground
  first_base@ps_hits_pos[[1]] <- 1L
  expect_identical(
    unname(ps_hits_pos(first_base, pos_shift = -200L)[[1]]),
    -200L
  )
})

test_that("direct scans populate complete foreground hit vectors", {
  motifs <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
  promoters <- readRDS(system.file("extdata", "prom_seq.rds", package = "PscanR"))[1:3]

  scanned <- ps_scan(as(motifs[[1]], "PSMatrix"), promoters)

  expect_identical(ps_fg_size(scanned), 3L)
  expect_identical(length(ps_hits_oligo(scanned)), 3L)
  expect_true(methods::validObject(scanned, test = TRUE))
})

test_that("results table labels and Z-score formula match implementation", {
  foreground <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  results <- PSMatrixList(foreground)
  names(results) <- TFBSTools::ID(foreground)
  table <- ps_results_table(results)
  scores <- ps_hits_score(foreground)
  n <- sum(!is.na(scores))
  expected_z <- (mean(scores, na.rm = TRUE) - ps_bg_avg(foreground)) /
    (ps_bg_std_dev(foreground) / sqrt(n))

  expect_identical(table$NAME[[1]], TFBSTools::name(foreground))
  expect_identical(row.names(table), TFBSTools::ID(foreground))
  expect_equal(table$ZSCORE[[1]], expected_z)
})
