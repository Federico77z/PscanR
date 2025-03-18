txdb <- txdbmaker::makeTxDbFromUCSC(genome="dm6", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:7] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 950, downstream = 50, use.names = TRUE)
prom_rng <- GenomicRanges::trim(prom_rng) 
prom_seq <- Biostrings::getSeq(x = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, prom_rng) #promoter sequences

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "insects"

J2022 <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts) #core Jaspar 2022 profiles for insects

J2022_PSBG <- PscanR::ps_build_bg(prom_seq, J2022, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2020_PSBG, "J2022_dm6_950u_50d_UCSC.psbg.txt")