TAIR_path <- system.file('extdata/BG_scripts', 'TAIR9_GFF3_genes.gff', package = 'PscanR')
txdb <- txdbmaker::makeTxDbFromGFF(TAIR_path, format = "gff3", dataSource="TAIR9", organism="Arabidopsis thaliana")
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:5]

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 450, downstream = 50, use.names = TRUE) 
GenomeInfoDb::seqlengths(prom_rng) <- GenomeInfoDb::seqlengths(BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9)[1:5]
prom_rng <- GenomicRanges::trim(prom_rng)
prom_seq <- Biostrings::getSeq(x = BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9, 
                               prom_rng)
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "plants"

J2020 <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts) #core Jaspar 2020 profiles for plants

J2020_PSBG <- PscanR::ps_build_bg(prom_seq, J2020, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2020_PSBG, "J2020_TAIR9_450u_50d_TAIR.psbg.txt")