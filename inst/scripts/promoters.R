# The aim of this script is to explain how the promoters.rds file was obtained.
# This file contains the promoter regions information for Homo sapiens (gene 
# assembly hg38), retrieved from the UCSC database. 
# The investigated promoter regions are 200bp downstream and 50bp upstream
# in respect to the trasncription start site. 

# Retrieve gene information from UCSC
txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated")
# Use only canonical chromosomes 
seqlevels(txdb) <- seqlevels(txdb)[1:24] 
# Compute the promoter regions 
prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE)
