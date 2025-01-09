## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
opts <- list()
opts[["collection"]] <- "CORE" # Core motif collection
opts[["tax_group"]] <- "vertebrates" # Vertebrate motifs
J2020 <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) # Download motif

file_path <- system.file("extdata", "example_file_nfkb100.txt", 
                         package = "PscanR")
BG_path <- system.file("extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt", 
                       package = "PscanR")
# Get the target regions
target <- read.csv(file_path, header = FALSE)

#txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", 
#                                    tablename="ncbiRefSeqCurated")

txdb_path <- system.file("extdata", "txdb.db", package = "PscanR")
txdb <- AnnotationDbi::loadDb(txdb_path)
seqlevels(txdb) <- seqlevels(txdb)[1:24]

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE)
prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)
target_prom_rng <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]

# Get the target sequences 
prom_seq <- Biostrings::getSeq(
  x = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
  target_prom_rng) 

# Build the background
J2020_PSBG <- PscanR::ps_build_bg_from_file(BG_path, J2020)

# Run PscanR
results <- PscanR::pscan(prom_seq, J2020_PSBG, 
                 BPPARAM = BiocParallel::MulticoreParam(1)) 
# Use SnowParam on Windows. 

# Create the table of results
table <- PscanR::ps_results_table(results) 
head(table)

