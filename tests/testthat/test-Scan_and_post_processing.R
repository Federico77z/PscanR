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
