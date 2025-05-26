txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") # Import GTF annotation from UCSC for Human genome 
# You can specify the organism of interest with the 'genome' parameter. 
# 'hg38' for Human, 'mm10' for Mouse, ... 
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:24] # Use only annotation for canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)

file_path <- system.file("extdata", "nrf1100.txt", package = "PscanR")
target <- read.csv(file_path, header = FALSE) 

target_prom_rng <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]
prom_seq <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                               target_prom_rng) 

saveRDS(prom_seq, "inst/extdata/prom_seq.rds")
