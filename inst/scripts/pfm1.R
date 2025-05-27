# The aim of this script is to explain how the pfm1.RData was generate. pfm1 
# corresponds to the first matrix of the PSMatrixList obtained as a result of 
# the PscanR algorithm. Pscan scans a DNAStringSet object with the Position
# Weight Matrices retrieved from the JASPAR database and adjusted based on the 
# background. 

# Load the pre-computed promoter sequences
file_path <- system.file("extdata", "prom_seq.rds", package = "PscanR")
prom_seq <- readRDS(file_path)
# We take only 50 sequences to reduce the computation
# time and the dimension of the final dataset. 
 
# Generate the background PSMatrixList
J2020_PSBG <- generate_psmatrixlist_from_background('JASPAR2020', 'hs', c(-200,+50), 'hg38')

# Execute the Pscan algorithm and view the result table

results <- pscan(prom_seq, J2020_PSBG, 
                BPPARAM = BiocParallel::SnowParam(1))
# Use MulticoreParam() for Unix systems (See BiocParallel package).
 
pfm1 <- results[['MA0506.1']]
saveRDS(pfm1, 'inst/extdata/pfm1.rds')
