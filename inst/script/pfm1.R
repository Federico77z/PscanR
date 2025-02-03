file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
prom_seq <- readRDS(file_path)
prom_seq <- prom_seq[1:50]
 
# Load JASPAR motif matrices for vertebrates
J2020_path <- system.file("extdata", "J2020.rda", package = "PscanR")
load(J2020_path)
 
bg_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
                        package = "PscanR")
J2020_PSBG <- ps_retrieve_bg_from_file(bg_path, J2020)

# Execute the Pscan algorithm and view the result table

results <- pscan(prom_seq, J2020_PSBG, 
                BPPARAM = BiocParallel::SnowParam(1))
# Use MulticoreParam() for Unix systems (See BiocParallel package).
 
pfm1 <- results[[1]]
save(pfm1, file = 'pfm1.RData')
