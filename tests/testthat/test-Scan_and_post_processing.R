test_that("PscanR works", {
  x <- Biostrings::DNAStringSet(c("ATGCTGCAATCGA", "CATGCTAAGCTAT", "GTACTACTAAATG",
                                  "TCAGACCATTAAA"))
  names(x) <- c("NM_001078.4","NM_000639.3", "NM_000756.8", "NM_001094.2")
  PFM1 <- PFMatrix(ID = "PSM1", name = "Example1", matrixClass = "PWM", 
                   profileMatrix = matrix(c(4, 19, 0, 0, 0, 0,
                                  16, 0, 20, 0, 0, 0, 
                                  0, 1, 0, 20, 0, 20,
                                  0, 0, 0, 0, 20, 0), 
                                nrow = 4, byrow = TRUE,
                                dimnames = list(c("A", "C", "G", "T"))))
  PFM2 <- PFMatrix(ID = "PSM2", name = "Example2", matrixClass = "PWM",
                   profileMatrix = matrix(c(3, 0, 0, 0, 0, 0,
                                  8, 0, 23, 0, 0, 0, 
                                  2, 23, 0, 23, 0, 24,
                                  11, 1, 1, 1, 24, 0), 
                                nrow = 4, byrow = TRUE,
                                dimnames = list(c("A", "C", "G", "T"))))
  PSM1 <- PSMatrix(PFM1, ps_bg_avg = 0.8267, ps_fg_avg = 0.8155478, 
                   ps_bg_std_dev = 0.07493493, ps_bg_size = 250L, 
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  PSM2 <- PSMatrix(PFM2, ps_bg_avg = 0.8806266, ps_fg_avg = 0.8679936, 
                   ps_bg_std_dev = 0.07161552, ps_bg_size = 250L,
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  pfms <- PSMatrixList(PSM1, PSM2)
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SerialParam())
  
  # The result is of PSMatrixList class
  
  expect_s4_class(result, "PFMatrixList")
  
  # The result is not empty 
  
  expect_gt(length(result), 0)
  
  # Invalid input 
  err_pfms <- c(5,6,89,4)
  expect_error(pscan(x, err_pfms, BPPARAM = BiocParallel::SerialParam()),
               "pfms is not an object of PSMatrixList class")
  
  
})

test_that("ps_result_table works", {
  x <- Biostrings::DNAStringSet(c("ATGCTGCAATCGA", "CATGCTAAGCTAT", "GTACTACTAAATG",
                                  "TCAGACCATTAAA"))
  names(x) <- c("NM_001078.4","NM_000639.3", "NM_000756.8", "NM_001094.2")
  PFM1 <- PFMatrix(ID = "PSM1", name = "Example1", matrixClass = "PWM", 
                   profileMatrix = matrix(c(4, 19, 0, 0, 0, 0,
                                            16, 0, 20, 0, 0, 0, 
                                            0, 1, 0, 20, 0, 20,
                                            0, 0, 0, 0, 20, 0), 
                                          nrow = 4, byrow = TRUE,
                                          dimnames = list(c("A", "C", "G", "T"))))
  PFM2 <- PFMatrix(ID = "PSM2", name = "Example2", matrixClass = "PWM",
                   profileMatrix = matrix(c(3, 0, 0, 0, 0, 0,
                                            8, 0, 23, 0, 0, 0, 
                                            2, 23, 0, 23, 0, 24,
                                            11, 1, 1, 1, 24, 0), 
                                          nrow = 4, byrow = TRUE,
                                          dimnames = list(c("A", "C", "G", "T"))))
  PSM1 <- PSMatrix(PFM1, ps_bg_avg = 0.8267, ps_fg_avg = 0.8155478, 
                   ps_bg_std_dev = 0.07493493, ps_bg_size = 250L, 
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  PSM2 <- PSMatrix(PFM2, ps_bg_avg = 0.8806266, ps_fg_avg = 0.8679936, 
                   ps_bg_std_dev = 0.07161552, ps_bg_size = 250L,
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  pfms <- PSMatrixList(PSM1, PSM2)
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SerialParam())
  
  table <- ps_results_table(result)
  filtered_table <- ps_results_table(result, FDR = max(table$FDR))
  empty_table <- ps_results_table(result, FDR = 0)
  
  # The type of the result
  
  expect_type(table, "list")
  expect_s3_class(filtered_table, "data.frame")
  expect_s3_class(empty_table, "data.frame")
  
  # The table is not empty 
  
  expect_gt(length(table), 0)
  expect_equal(nrow(filtered_table), nrow(table))
  expect_equal(nrow(empty_table), 0)
  
  # Invalid input 
  err_pfms <- c(5,6,89,4)
  expect_error(ps_results_table(err_pfms), 
               "pfms is not an object of PSMatrixList class")
  expect_error(ps_results_table(result, FDR = -0.1),
               "FDR must be a single numeric value between 0 and 1")
  expect_error(ps_results_table(result, FDR = 1.1),
               "FDR must be a single numeric value between 0 and 1")
  expect_error(ps_results_table(result, FDR = c(0.01, 0.05)),
               "FDR must be a single numeric value between 0 and 1")
  expect_error(ps_results_table(result, FDR = NA_real_),
               "FDR must be a single numeric value between 0 and 1")
  
})

test_that("ps_z_table works",{
  x <- Biostrings::DNAStringSet(c("ATGCTGCAATCGA", "CATGCTAAGCTAT", "GTACTACTAAATG",
                                  "TCAGACCATTAAA"))
  names(x) <- c("NM_001078.4","NM_000639.3", "NM_000756.8", "NM_001094.2")
  PFM1 <- PFMatrix(ID = "PSM1", name = "Example1", matrixClass = "PWM", 
                   profileMatrix = matrix(c(4, 19, 0, 0, 0, 0,
                                            16, 0, 20, 0, 0, 0, 
                                            0, 1, 0, 20, 0, 20,
                                            0, 0, 0, 0, 20, 0), 
                                          nrow = 4, byrow = TRUE,
                                          dimnames = list(c("A", "C", "G", "T"))))
  PFM2 <- PFMatrix(ID = "PSM2", name = "Example2", matrixClass = "PWM",
                   profileMatrix = matrix(c(3, 0, 0, 0, 0, 0,
                                            8, 0, 23, 0, 0, 0, 
                                            2, 23, 0, 23, 0, 24,
                                            11, 1, 1, 1, 24, 0), 
                                          nrow = 4, byrow = TRUE,
                                          dimnames = list(c("A", "C", "G", "T"))))
  PSM1 <- PSMatrix(PFM1, ps_bg_avg = 0.8267, ps_fg_avg = 0.8155478, 
                   ps_bg_std_dev = 0.07493493, ps_bg_size = 250L, 
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  PSM2 <- PSMatrix(PFM2, ps_bg_avg = 0.8806266, ps_fg_avg = 0.8679936, 
                   ps_bg_std_dev = 0.07161552, ps_bg_size = 250L,
                   ps_seq_names = c("NM_001078.4","NM_000639.3",
                                    "NM_000756.8", "NM_001094.2"))
  pfms <- PSMatrixList(PSM1, PSM2)
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SerialParam())
  z_table <- ps_z_table(result)
  
  # Expect a valid output
  
  expect_gt(length(z_table), 0)
  
  # The type of the result 
  
  expect_type(table, "closure")
  
  # Invalid input
  err_pfms <- c(5,6,89,4)
  expect_error(ps_z_table(err_pfms),
               "pfms is not an object of PSMatrixList class")

})

test_that(".ps_scan_s scans correctly and handles edge cases", {
  # Motif with the unambiguous consensus "ACGTAC" (reverse complement "GTACGT")
  PFM <- PFMatrix(ID = "CONS", name = "Consensus", matrixClass = "PWM",
                  profileMatrix = matrix(c(20,  0,  0,  0, 20,  0,
                                            0, 20,  0,  0,  0, 20,
                                            0,  0, 20,  0,  0,  0,
                                            0,  0,  0, 20,  0,  0),
                                         nrow = 4, byrow = TRUE,
                                         dimnames = list(c("A", "C", "G", "T"))))
  PSM <- PSMatrix(PFM, ps_bg_avg = 0.5, ps_fg_avg = 0.5,
                  ps_bg_std_dev = 0.1, ps_bg_size = 100L)

  W <- ncol(TFBSTools::Matrix(PSM))
  nr <- length(PscanR:::.PS_ALPHABET(PSM))
  M <- matrix(as.numeric(TFBSTools::Matrix(PSM)), nrow = nr, ncol = W)
  M_rc <- matrix(as.numeric(TFBSTools::Matrix(reverseComplement(PSM))),
                 nrow = nr, ncol = W)
  # Raw score of a perfect consensus match equals the matrix maximum score.
  max_raw <- as.numeric(Biostrings::maxScore(TFBSTools::Matrix(PSM)))

  # Forward consensus embedded at position 5
  fwd <- PscanR:::.ps_scan_s(PSM, "TTTTACGTACTTTT", M, M_rc, W)
  expect_equal(fwd$strand, "+")
  expect_equal(fwd$pos, 5L)
  expect_equal(fwd$oligo, "ACGTAC")
  expect_equal(fwd$score, max_raw)

  # Reverse-complement consensus embedded at position 5
  rev <- PscanR:::.ps_scan_s(PSM, "TTTTGTACGTTTTT", M, M_rc, W)
  expect_equal(rev$strand, "-")
  expect_equal(rev$pos, 5L)
  expect_equal(rev$oligo, "GTACGT")
  expect_equal(rev$score, max_raw)

  # Soft-masked (lower-case) bases are scored as their upper-case equivalent
  lc <- PscanR:::.ps_scan_s(PSM, "ttttacgtactttt", M, M_rc, W)
  expect_equal(lc$pos, 5L)
  expect_equal(lc$score, max_raw)

  # Windows overlapping an N are skipped; the clean consensus still wins
  ncase <- PscanR:::.ps_scan_s(PSM, "NNACGTACNN", M, M_rc, W)
  expect_equal(ncase$pos, 3L)
  expect_equal(ncase$oligo, "ACGTAC")
  expect_equal(ncase$score, max_raw)

  # Sequence shorter than the motif -> NA score (not -Inf), no error
  short <- PscanR:::.ps_scan_s(PSM, "ACG", M, M_rc, W)
  expect_true(is.na(short$score))

  # Fully N-masked sequence -> NA score, no error
  alln <- PscanR:::.ps_scan_s(PSM, "NNNNNNNN", M, M_rc, W)
  expect_true(is.na(alln$score))
})

test_that(".ps_scan_batched matches the per-sequence kernel bitwise", {
  PFM <- PFMatrix(ID = "CONS", name = "Consensus", matrixClass = "PWM",
                  profileMatrix = matrix(c(20,  0,  0,  0, 20,  0,
                                            0, 20,  0,  0,  0, 20,
                                            0,  0, 20,  0,  0,  0,
                                            0,  0,  0, 20,  0,  0),
                                         nrow = 4, byrow = TRUE,
                                         dimnames = list(c("A", "C", "G", "T"))))
  PSM <- PSMatrix(PFM, ps_bg_avg = 0.5, ps_fg_avg = 0.5,
                  ps_bg_std_dev = 0.1, ps_bg_size = 100L)
  W <- ncol(TFBSTools::Matrix(PSM))
  nr <- length(PscanR:::.PS_ALPHABET(PSM))
  M <- matrix(as.numeric(TFBSTools::Matrix(PSM)), nrow = nr, ncol = W)
  M_rc <- matrix(as.numeric(TFBSTools::Matrix(reverseComplement(PSM))),
                 nrow = nr, ncol = W)

  # Equal-length sequences covering: forward consensus, reverse-complement
  # consensus, soft-masked (lower-case), an internal N, and an all-N sequence.
  seqs <- c(fwd  = "TTTTACGTACTTTT",
            rev  = "TTTTGTACGTTTTT",
            lc   = "ttttacgtactttt",
            nmix = "TTTTACGTACTTNN",
            alln = "NNNNNNNNNNNNNN")
  S <- PscanR:::.ps_encode_seqs(unname(seqs))
  expect_false(is.null(S))                 # all equal length -> fast path
  expect_identical(dim(S), c(length(seqs), unname(nchar(seqs[1L]))))
  expect_true(attr(S, "has_ambiguous"))
  batched <- PscanR:::.ps_scan_batched(unname(seqs), S, M, M_rc, W)

  # Reference: the per-sequence kernel, one sequence at a time.
  ref <- lapply(unname(seqs), function(s) PscanR:::.ps_scan_s(PSM, s, M, M_rc, W))
  ref_score  <- vapply(ref, function(r) as.numeric(r$score),  numeric(1))
  ref_pos    <- vapply(ref, function(r) as.integer(r$pos),    integer(1))
  ref_strand <- vapply(ref, function(r) as.character(r$strand), character(1))
  ref_oligo  <- vapply(ref, function(r) as.character(r$oligo), character(1))

  expect_identical(batched$score,  ref_score)   # bitwise
  expect_identical(batched$pos,    ref_pos)
  expect_identical(batched$strand, ref_strand)
  expect_identical(batched$oligo,  ref_oligo)

  # The no-ambiguity fast path skips NA-mask allocation but remains bitwise
  # identical to the per-sequence kernel.
  plain_seqs <- unname(seqs[c("fwd", "rev", "lc")])
  plain_encoded <- PscanR:::.ps_encode_seqs(plain_seqs)
  expect_false(attr(plain_encoded, "has_ambiguous"))
  expect_s4_class(attr(plain_encoded, "subject"), "DNAString")
  plain_batched <- PscanR:::.ps_scan_batched(
    plain_seqs, plain_encoded, M, M_rc, W
  )
  plain_ref <- lapply(plain_seqs, function(s) {
    PscanR:::.ps_scan_s(PSM, s, M, M_rc, W)
  })
  expect_identical(
    plain_batched$score,
    vapply(plain_ref, function(r) as.numeric(r$score), numeric(1))
  )
  expect_identical(
    plain_batched$pos,
    vapply(plain_ref, function(r) as.integer(r$pos), integer(1))
  )
  expect_identical(
    plain_batched$strand,
    vapply(plain_ref, function(r) as.character(r$strand), character(1))
  )
  expect_identical(
    plain_batched$oligo,
    vapply(plain_ref, function(r) as.character(r$oligo), character(1))
  )

  # Equal-length sequences all shorter than the motif -> NA, matching the kernel.
  shortseqs <- c(a = "ACG", b = "TGC")
  Ss <- PscanR:::.ps_encode_seqs(unname(shortseqs))
  bshort <- PscanR:::.ps_scan_batched(unname(shortseqs), Ss, M, M_rc, W)
  expect_true(all(is.na(bshort$score)))
})

test_that("p-values use the exact upper tail and do not underflow", {
  x <- Biostrings::DNAStringSet(c("ATGCTGCAATCGA", "CATGCTAAGCTAT",
                                  "GTACTACTAAATG", "TCAGACCATTAAA"))
  names(x) <- c("NM_001078.4", "NM_000639.3", "NM_000756.8", "NM_001094.2")
  PFM1 <- PFMatrix(ID = "PSM1", name = "Example1", matrixClass = "PWM",
                   profileMatrix = matrix(c(4, 19, 0, 0, 0, 0,
                                            16, 0, 20, 0, 0, 0,
                                            0, 1, 0, 20, 0, 20,
                                            0, 0, 0, 0, 20, 0),
                                          nrow = 4, byrow = TRUE,
                                          dimnames = list(c("A", "C", "G", "T"))))
  probe <- PSMatrix(PFM1, ps_bg_avg = 0.5, ps_fg_avg = NA_real_,
                    ps_bg_std_dev = 0.05, ps_bg_size = 250L,
                    ps_seq_names = names(x))
  probed <- pscan(x, PSMatrixList(probe), BPPARAM = BiocParallel::SerialParam())
  fg <- ps_fg_avg(probed[[1]])
  n <- length(ps_hits_score(probed[[1]]))

  # Choose a background so the z-score lands past the point where the old
  # 1 - pnorm(z) formulation collapses to exactly zero (z > 8.3), but well
  # short of where the true upper tail itself underflows (z ~ 37).
  target_z <- 12
  bg_avg <- fg - 0.05
  bg_sd <- (fg - bg_avg) * sqrt(n) / target_z
  enriched <- PSMatrix(PFM1, ps_bg_avg = bg_avg, ps_fg_avg = NA_real_,
                       ps_bg_std_dev = bg_sd, ps_bg_size = 250L,
                       ps_seq_names = names(x))
  result <- pscan(x, PSMatrixList(enriched),
                  BPPARAM = BiocParallel::SerialParam())
  z <- ps_zscore(result[[1]])
  p <- ps_pvalue(result[[1]])

  expect_equal(unname(z), target_z, tolerance = 1e-8)
  # The bug this guards against: 1 - pnorm(z) is exactly 0 here.
  expect_identical(1 - stats::pnorm(unname(z)), 0)
  expect_gt(p, 0)
  expect_identical(p, stats::pnorm(unname(z), lower.tail = FALSE))

  # The z statistic keeps its historical "z" name.
  expect_identical(names(z), "z")

  # The invariant holds for an ordinary, non-extreme motif too.
  ordinary <- ps_pvalue(probed[[1]])
  expect_identical(
    ordinary,
    stats::pnorm(unname(ps_zscore(probed[[1]])), lower.tail = FALSE)
  )
})

test_that("pscan_fullBG resolves the requested splice variant, not its gene", {
  # Two splice variants of one Arabidopsis gene with *different* promoter
  # sequences. Before transcript resolution became scheme-aware, both
  # identifiers were stripped to "AT1G01110" and match() returned whichever
  # came first, so the caller could silently receive the other variant's hits.
  variants <- Biostrings::DNAStringSet(c(
    AT1G01110.1 = "ATGCTGCAATCGATTTTTTTTTT",
    AT1G01110.2 = "GGGGGGGGGGCATGCTAAGCTAT",
    AT1G01160.1 = "GTACTACTAAATGCCCCCCCCCC",
    AT1G01200.1 = "TCAGACCATTAAAGGGGGGGGGG"
  ))
  PFM1 <- PFMatrix(
    ID = "PSM1", name = "Example1", matrixClass = "PWM",
    profileMatrix = matrix(
      c(4, 19, 0, 0, 0, 0,
        16, 0, 20, 0, 0, 0,
        0, 1, 0, 20, 0, 20,
        0, 0, 0, 0, 20, 0),
      nrow = 4, byrow = TRUE, dimnames = list(c("A", "C", "G", "T"))
    )
  )
  pfms <- PFMatrixList(PFM1)

  full_bg <- ps_build_bg(
    variants, pfms, BPPARAM = BiocParallel::SerialParam(), fullBG = TRUE
  )

  bg_scores <- ps_hits_score_bg(full_bg[[1]])
  # The test is only meaningful if the two variants score differently.
  expect_false(identical(
    unname(bg_scores[["AT1G01110.1"]]), unname(bg_scores[["AT1G01110.2"]])
  ))

  # Ask for .2 and deliberately not .1. (Three identifiers are the minimum the
  # z-test accepts.)
  wanted <- c("AT1G01110.2", "AT1G01160.1", "AT1G01200.1")
  retrieved <- pscan_fullBG(wanted, full_bg, quiet = TRUE)

  expect_identical(ps_seq_names(retrieved[[1]]), wanted)
  expect_identical(
    unname(ps_hits_score(retrieved[[1]])), unname(bg_scores[wanted])
  )
  # Specifically: the score returned for .2 is .2's, not .1's.
  expect_identical(
    unname(ps_hits_score(retrieved[[1]])[[1L]]),
    unname(bg_scores[["AT1G01110.2"]])
  )

  # The bare gene identifier is not a transcript here, so it must be reported
  # as absent rather than silently resolving to one of its variants.
  expect_warning(
    gene_query <- pscan_fullBG(
      c("AT1G01110", wanted), full_bg, quiet = TRUE
    ),
    "absent from the background"
  )
  expect_identical(ps_seq_names(gene_query[[1]]), wanted)
})

# ---------------------------------------------------------------------------
# Motif class and the ranked bar plot
# ---------------------------------------------------------------------------

scan_bundled_motifs <- function() {
  prom_seq <- readRDS(
    system.file("extdata", "prom_seq.rds", package = "PscanR")
  )[1:10]
  J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
  bg <- ps_retrieve_bg_from_file(
    system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
                package = "PscanR"),
    J2020
  )
  bg <- bg[c("MA0506.1", "MA0632.2", "MA0611.1",
             "MA0685.1", "MA0698.1", "MA0699.1")]
  pscan(prom_seq, bg, BPPARAM = BiocParallel::SerialParam())
}

test_that("ps_motif_class reads the class off the matrices", {
  results <- scan_bundled_motifs()
  classes <- ps_motif_class(results)

  expect_length(classes, length(results))
  # ID() is itself a named vector, so compare against its values.
  expect_identical(names(classes), unname(TFBSTools::ID(results)))
  expect_type(classes, "character")
  expect_identical(
    unname(classes[c("MA0506.1", "MA0632.2")]),
    c("Basic leucine zipper factors (bZIP)",
      "Basic helix-loop-helix factors (bHLH)")
  )
  # Every motif gets exactly one label, so it is safe as a grouping variable.
  expect_false(any(is.na(classes)))
  expect_true(all(nzchar(classes)))

  expect_error(ps_motif_class(data.frame(a = 1)), "PSMatrixList")
})

test_that("ps_motif_class reports missing and multiple classes usefully", {
  J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
  one <- J2020[["MA0506.1"]]

  # The JASPAR sentinel for "no class" is the empty string, not NA.
  blank <- one
  blank@matrixClass <- ""
  expect_identical(unname(ps_motif_class(PFMatrixList(blank))), "Unclassified")

  # A few JASPAR matrices carry two classes; only the first is reported so the
  # result stays one label per motif.
  paired <- one
  paired@matrixClass <- c("First class", "Second class")
  expect_identical(unname(ps_motif_class(PFMatrixList(paired))), "First class")
})

test_that("ps_motif_class strips the stray whitespace JASPAR ships", {
  J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
  one <- J2020[["MA0506.1"]]

  # JASPAR 2020 gives IRF3, IRF4 and IRF5 a trailing space, which used to split
  # "Tryptophan cluster factors" into two grouping levels that plot as two
  # legend entries with the same label.
  clean <- one
  clean@matrixClass <- "Tryptophan cluster factors"
  padded <- one
  padded@matrixClass <- "Tryptophan cluster factors "

  expect_identical(
    unname(ps_motif_class(PFMatrixList(clean, padded))),
    rep("Tryptophan cluster factors", 2L)
  )

  # Whitespace-only is as classless as the empty string.
  whitespace <- one
  whitespace@matrixClass <- "   "
  expect_identical(
    unname(ps_motif_class(PFMatrixList(whitespace))), "Unclassified"
  )
})

test_that("ps_motif_barplot returns a composable ggplot", {
  results <- scan_bundled_motifs()
  plot <- ps_motif_barplot(results, n = 4)

  expect_s3_class(plot, "ggplot")
  expect_identical(nrow(plot$data), 4L)
  # Nothing is drawn until printing, so further layers can still be added.
  expect_s3_class(plot + ggplot2::labs(title = "t"), "ggplot")
})

test_that("ps_motif_barplot ranks according to the statistic", {
  results <- scan_bundled_motifs()
  res_table <- ps_results_table(results)

  by_z <- ps_motif_barplot(results, n = 6, statistic = "ZSCORE")
  expect_identical(
    as.character(by_z$data$motif[1]),
    res_table$NAME[which.max(res_table$ZSCORE)]
  )

  # P-values rank the other way and are plotted on a log scale, because a bar
  # chart of raw p-values spanning many orders of magnitude is unreadable.
  by_p <- ps_motif_barplot(results, n = 6, statistic = "P.VALUE")
  expect_identical(
    as.character(by_p$data$motif[1]),
    res_table$NAME[which.min(res_table$P.VALUE)]
  )
  expect_equal(by_p$data$value[1], -log10(min(res_table$P.VALUE)))
  expect_identical(by_p$labels$x, "-log10(P.VALUE)")

  expect_error(ps_motif_barplot(results, statistic = "NOPE"), "arg")
})

test_that("ps_motif_barplot groups by class or by a supplied vector", {
  results <- scan_bundled_motifs()

  by_class <- ps_motif_barplot(results, n = 6, group = "class")
  expect_identical(
    sort(unique(as.character(by_class$data$group))),
    sort(unique(unname(ps_motif_class(results))))
  )

  supplied <- rep(c("a", "b"), length.out = length(results))
  by_vector <- ps_motif_barplot(results, n = 6, group = supplied)
  expect_setequal(unique(as.character(by_vector$data$group)), c("a", "b"))

  # A named vector is matched by identifier rather than by position.
  named <- stats::setNames(supplied, TFBSTools::ID(results))
  expect_s3_class(ps_motif_barplot(results, n = 6, group = named), "ggplot")

  expect_error(
    ps_motif_barplot(results, group = c("only", "two")),
    "one entry per motif"
  )
  expect_error(ps_motif_barplot(results, group = 1:6), "character or")
})

test_that("ps_motif_barplot clamps n and honours the FDR filter", {
  results <- scan_bundled_motifs()

  expect_identical(
    nrow(ps_motif_barplot(results, n = 99)$data), length(results)
  )
  expect_error(ps_motif_barplot(results, n = 0), "positive")

  res_table <- ps_results_table(results)
  cutoff <- stats::median(res_table$FDR)
  filtered <- ps_motif_barplot(results, n = 99, FDR = cutoff)
  expect_identical(nrow(filtered$data), sum(res_table$FDR <= cutoff))
  expect_error(ps_motif_barplot(results, FDR = 0), "No motifs left")
})

test_that("ps_motif_barplot accepts a results table when told the grouping", {
  results <- scan_bundled_motifs()
  res_table <- ps_results_table(results)

  expect_s3_class(ps_motif_barplot(res_table, n = 4), "ggplot")
  # A table carries no motif metadata, so "class" cannot be derived from it.
  expect_error(
    ps_motif_barplot(res_table, group = "class"), "results table does not"
  )
  expect_error(ps_motif_barplot(list(1, 2)), "PSMatrixList or a data.frame")
})

# ---------------------------------------------------------------------------
# Density plots
# ---------------------------------------------------------------------------

test_that("ps_density_plot returns a composable ggplot over the density grid", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  plot <- ps_density_plot(pfm, shift = -200)

  expect_s3_class(plot, "ggplot")
  # The curve is the grid of the density object, not a recomputed geom_density.
  positions <- ps_hits_pos(pfm, pos_shift = -200)[
    ps_hits_score(pfm) >= ps_bg_avg(pfm)
  ]
  expected <- stats::density(
    positions,
    from = min(positions), to = max(positions)
  )
  expect_identical(nrow(plot$data), length(expected$x))
  expect_equal(plot$data$x, expected$x)
  expect_equal(plot$data$y, expected$y)

  # Nothing is drawn until printing, so further layers can still be added.
  expect_s3_class(plot + ggplot2::labs(title = "t"), "ggplot")
})

test_that("ps_density_plot does not draw beyond the promoter", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  plot <- ps_density_plot(pfm, shift = -200)
  positions <- ps_hits_pos(pfm, pos_shift = -200)[
    ps_hits_score(pfm) >= ps_bg_avg(pfm)
  ]

  # stats::density() evaluates three bandwidths past the data by default, which
  # would draw hits outside the window that was scanned.
  expect_equal(min(plot$data$x), min(positions))
  expect_equal(max(plot$data$x), max(positions))

  # ... and the default would have: it runs past the data at both ends.
  unbounded <- stats::density(positions)
  expect_lt(min(unbounded$x), min(positions))
  expect_gt(max(unbounded$x), max(positions))
})

test_that("ps_density_plot corrects the estimate when given a window", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  win <- c(-200, 50)

  clipped <- ps_density_plot(pfm, shift = -200)
  reflected <- ps_density_plot(pfm, shift = -200, window = win)

  # The corrected curve spans every position a hit could be reported at, which
  # stops one motif width short of the end of the window.
  expect_equal(
    range(reflected$data$x), c(win[[1]] + 1, win[[2]] - ncol(pfm) + 1)
  )

  # Reflection returns the mass that a naive estimate lets escape, so it
  # integrates to one over the window while the clipped one falls short.
  area <- function(p) sum(diff(p$data$x)[[1]] * p$data$y)
  expect_lt(area(clipped), 0.99)
  expect_equal(area(reflected), 1, tolerance = 0.01)

  # ... and it is higher at the edges, which is where the mass went.
  expect_gt(reflected$data$y[[1]], clipped$data$y[[1]])
})

test_that("ps_density_plot validates the window", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))

  expect_error(
    ps_density_plot(pfm, shift = -200, window = c(-100, 0)),
    "fall outside"
  )
  expect_error(
    ps_density_plot(pfm, shift = -200, window = c(0, 0)),
    "two finite, distinct values"
  )
  # Narrower than the motif, so no hit could be reported in it at all.
  expect_error(
    ps_density_plot(pfm, shift = -200, window = c(-200, -197)),
    "too short to report a hit"
  )
})

test_that("ps_density_plot shifts positions relative to the TSS", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))

  unshifted <- ps_density_plot(pfm)
  shifted <- ps_density_plot(pfm, shift = -200)

  expect_equal(shifted$data$x, unshifted$data$x - 200)
})

test_that("ps_density_plot honours the score threshold", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))

  # The promoter count is reported in the title, so the threshold is visible
  # in the returned object without having to render it.
  promoters <- function(plot) {
    as.integer(sub("^.* across ([0-9]+) promoter.*$", "\\1", plot$labels$title))
  }
  all_hits <- promoters(ps_density_plot(pfm, st = "all"))
  loose <- promoters(ps_density_plot(pfm, st = "loose"))
  strict <- promoters(ps_density_plot(pfm, st = "strict"))

  expect_identical(all_hits, length(ps_hits_score(pfm)))
  expect_gt(all_hits, loose)
  expect_gt(loose, strict)

  # An unrecognised threshold falls back to "loose" rather than failing.
  expect_warning(
    fallback <- ps_density_plot(pfm, st = "nonsense"),
    "reverting to loose"
  )
  expect_equal(fallback$data, ps_density_plot(pfm, st = "loose")$data)
})

test_that("ps_density_distances_plot returns a ggplot and checks its inputs", {
  results <- scan_bundled_motifs()
  plot <- ps_density_distances_plot(
    results[["MA0506.1"]], results[["MA0632.2"]], "all", "loose"
  )

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot + ggplot2::labs(title = "t"), "ggplot")
  expect_match(plot$labels$x, "Distances between")

  expect_error(
    ps_density_distances_plot(results[["MA0506.1"]], "not a matrix"),
    "Both object must be of class PSMatrix"
  )
})

test_that("ps_density_plot can draw a binned profile on the density scale", {
  pfm <- readRDS(system.file("extdata", "pfm1.rds", package = "PscanR"))
  plot <- ps_density_plot(pfm, shift = -200, window = c(-200, 50), bins = 12)

  bars <- ggplot2::layer_data(plot, 1)
  expect_identical(nrow(bars), 12L)

  # Bars and curve share one axis because both integrate to one -- no
  # secondary axis and no rescaling factor.
  expect_equal(sum(bars$y * (bars$xmax - bars$xmin)), 1, tolerance = 1e-6)
  expect_equal(
    sum(diff(plot$data$x)[[1]] * plot$data$y), 1,
    tolerance = 0.01
  )

  # The bars span the same support the curve does.
  expect_equal(min(bars$xmin), min(plot$data$x))
  expect_equal(max(bars$xmax), max(plot$data$x))

  # Without bins the first layer is the filled area, not bars.
  plain <- ps_density_plot(pfm, shift = -200, window = c(-200, 50))
  expect_lt(nrow(ggplot2::layer_data(plain, 1)), 512 + 12)
  expect_error(ps_density_plot(pfm, bins = 0), "positive")
})

test_that("the binned profile counts every value exactly once", {
  set.seed(1)
  x <- runif(500, -100, 0)
  bars <- PscanR:::.ps_binned_profile(x, bins = 10, limits = c(-100, 0))

  # Recovering counts from density heights must return the whole sample.
  counts <- bars$density * length(x) * bars$width
  expect_equal(sum(counts), length(x))
  expect_true(all(counts >= 0))
})

test_that("ps_hit_score_plot pairs each promoter's position with its score", {
  results <- scan_bundled_motifs()
  one <- results[["MA0506.1"]]
  plot <- ps_hit_score_plot(one, shift = -200)

  expect_s3_class(plot, "ggplot")
  # One point per promoter, carrying both coordinates.
  expect_identical(nrow(plot$data), length(ps_hits_score(one)))
  expect_equal(plot$data$score, unname(ps_hits_score(one)))
  expect_equal(plot$data$position, unname(ps_hits_pos(one, pos_shift = -200)))

  # The bands are cut by the two thresholds `st` names elsewhere, so the
  # counts must agree with what those thresholds keep.
  loose <- ps_bg_avg(one)
  strict <- loose + ps_bg_std_dev(one)
  counts <- table(plot$data$band)
  expect_identical(
    unname(counts[["above mean + SD (strict)"]]),
    sum(ps_hits_score(one) >= strict)
  )
  expect_identical(
    unname(counts[["above mean (loose)"]]),
    sum(ps_hits_score(one) >= loose & ps_hits_score(one) < strict)
  )
  expect_identical(
    unname(counts[["below background mean"]]),
    sum(ps_hits_score(one) < loose)
  )
  # Every promoter lands in exactly one band.
  expect_identical(sum(counts), nrow(plot$data))
  expect_false(any(is.na(plot$data$band)))
})

test_that("ps_hit_score_plot panels a PSMatrixList in the given order", {
  results <- scan_bundled_motifs()
  ids <- c("MA0632.2", "MA0506.1")
  pair <- results[match(ids, TFBSTools::ID(results))]
  pair <- BiocGenerics::do.call(PSMatrixList, as.list(pair))

  plot <- ps_hit_score_plot(pair, shift = -200)
  expect_identical(
    levels(plot$data$motif), unname(TFBSTools::name(pair))
  )
  expect_identical(
    nrow(plot$data),
    sum(vapply(as.list(pair), function(m) length(ps_hits_score(m)), integer(1)))
  )

  expect_error(ps_hit_score_plot(data.frame(a = 1)), "PSMatrix")
})
