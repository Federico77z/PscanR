# The aim of this script is to explain how the mega_pfm1.RData was generate. 
# mega_pfm1 corresponds to the first matrix of the background PSMatrixList 
# obtained as a result of the PscanR algorithm applied to all the promoters 
# regions in a specified organism, scanned with the PFM retrieved from the 
# JASPAR database
txdb <- txdbmaker::makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:24] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
prom_rng <- prom_rng[1:50]
prom_seq <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, prom_rng) #promoter sequences

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2020 <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) #core Jaspar 2020 profiles for vertebrates

J2020_PSBG <- PscanR::ps_build_bg(prom_seq, J2020, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

mega_pfm1 <- J2020_PSBG[[1]]

save(mega_pfm1, file = "mega_pfm1.RData")