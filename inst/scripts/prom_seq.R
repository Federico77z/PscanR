txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") # Import GTF annotation from UCSC for Human genome 
# You can specify the organism of interest with the 'genome' parameter. 
# 'hg38' for Human, 'mm10' for Mouse, ... 
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:24] # Use only annotation for canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE)

file_path <- system.file("extdata", "nrf1100.txt", package = "PscanR")
target <- read.csv(file_path, header = FALSE)

# nrf1100.txt carries RefSeq versions and is matched exactly, which pins the
# fixture to the annotation release it was built from. That is deliberate: a
# version bump can move the transcription start site, so matching on the bare
# accession would quietly return a different promoter sequence under the same
# name. Fewer than 90 promoters here means an accession has been versioned
# upstream and the target list needs revisiting, not that this script is wrong.
target_prom_rng <- prom_rng[prom_rng$tx_name %in% target[,1]]
prom_seq <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                               target_prom_rng) 

saveRDS(prom_seq, "inst/extdata/prom_seq.rds")
