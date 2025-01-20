txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") # Import GTF annotation from UCSC for Human genome 
# You can specify the organism of interest with the 'genome' parameter. 
# 'hg38' for Human, 'mm10' for Mouse, ... 
seqlevels(txdb) <- seqlevels(txdb)[1:24] # Use only annotation for canonical chromosomes

prom_rng <- promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
prom_seq <- Biostrings::getSeq(
  x = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
  prom_rng) 

saveRDS(prom_seq, "inst/extdata/prom_seq.rds")