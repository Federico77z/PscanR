txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")
seqlevels(txdb) <- seqlevels(txdb)[1:24]
prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE)
