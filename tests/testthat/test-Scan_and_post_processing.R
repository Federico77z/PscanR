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
  
  # The type of the result
  
  expect_type(table, "list")
  
  # The table is not empty 
  
  expect_gt(length(table), 0)
  
  # Invalid input 
  err_pfms <- c(5,6,89,4)
  expect_error(ps_results_table(err_pfms), 
               "pfms is not an object of PSMatrixList class")
  
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