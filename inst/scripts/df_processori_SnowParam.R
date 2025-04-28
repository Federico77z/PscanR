# This script explain how the df_processori.rds file was obtained. 
# This file contains a data frame with 8 columns, corresponding to the different
# number of cores. Rows represent the number of iteration. Each row contains the 
# evaluated time for each core using the SnowParam() function. 

# Import gtf annotation from UCSC
txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")

# Select only canonical chromosomes 
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:24] 

# Background generation. JASPAR core dataset for Vertebrates, 2020 release
J2020_PSBG <- PscanR::generate_psmatrixlist_from_background( 'Jaspar2020', 'hs', c(-950,50), 'hg38') 

# Retrieve promoter regions 
prom_rng <- GenomicFeatures::promoters(txdb, upstream = 950, downstream = 50, use.names = TRUE) 

# Version numbers removal from transcript IDs
prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name) 

# Retrieve target sequence identifiers. You can download it from the official 
# Pscan platform
file_path <- system.file('extdata', 'liver.txt', package = 'PscanR')
target <- read.csv(file_path, header = F) 

# Filters promoter regions corresponding to the target gene list
prom_range <- prom_rng[prom_rng$tx_name_clean %in% target[,1]] 

# Extract sequences
prom_seq <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38, prom_range)

iterations <- 100 # Number of iterations
core_numbers <- c(1,2,4,6,8,12,18,24) # Cores used

execution_time <- function(prom_seq, bg, ncores, output_file){
  temp <- numeric(iterations)
  
  for (i in 1:iterations){
    timing <- system.time({
      results <- PscanR::pscan(prom_seq, bg, BPPARAM = BiocParallel::SnowParam(ncores))}) 
    if (i == 1){
      saveRDS(PscanR::ps_results_table(results), output_file) 
      # Analyse all the files to be sure that results don't change when changing the core numbers
    }
    temp[i] <- timing["elapsed"]
  }
  return(temp)
}

all_cores <- function(core_numbers, prom_seq, bg, iterations) {
  results_list <- list()
  
  for (core in core_numbers) {
    output_file <- paste0("core", core, "_liver_950d_50u.rds")
    results_list[[as.character(core)]] <- execution_time(prom_seq, bg, core, output_file)
  }
  
  return(as.data.frame(results_list))
} 

df_processori <- all_cores(core_numbers, prom_seq, bg_hg38_2020, iterations)
colnames(df_processori) <- core_numbers

saveRDS(df_processori, "df_processori_SnowParam.rds")