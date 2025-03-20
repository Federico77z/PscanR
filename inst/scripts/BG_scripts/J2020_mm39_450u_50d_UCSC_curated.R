txdb <- txdbmaker::makeTxDbFromUCSC(genome="mm39", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:21] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 450, downstream = 50, use.names = TRUE) 
prom_seq <- Biostrings::getSeq(x = BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39, 
                               prom_rng) #promoter sequences
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2020 <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) #core Jaspar 2020 profiles for vertebrates

J2020_PSBG <- PscanR::ps_build_bg(prom_seq, J2020, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2020_PSBG, "J2020_mm39_450u_50d_UCSC.psbg.txt")