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
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SnowParam(1))
  
  # The result is of PSMatrixList class
  
  expect_s4_class(result, "PFMatrixList")
  
  # The result is not empty 
  
  expect_gt(length(result), 0)
  
  # Invalid input 
  err_pfms <- c(5,6,89,4)
  expect_error(pscan(x, err_pfms, BPPARAM = BiocParallel::SnowParam(1)), 
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
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SnowParam(1))
  
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
  
  result <- pscan(x, pfms, BPPARAM = BiocParallel::SnowParam(1))
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
