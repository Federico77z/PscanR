TAIR_path <- system.file('extdata/BG_files', 'TAIR9_GFF3_genes.gff', package = 'PscanR')
txdb <- txdbmaker::makeTxDbFromGFF(TAIR_path, format = "gff3", dataSource="TAIR9", organism="Arabidopsis thaliana")
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:5]

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
GenomeInfoDb::seqlengths(prom_rng) <- GenomeInfoDb::seqlengths(BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9)[1:5]
prom_rng <- GenomicRanges::trim(prom_rng)
prom_seq <- Biostrings::getSeq(x = BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9, 
                               prom_rng)
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "plants"

httr::set_config(httr::config(ssl_verifypeer = 0L))
JASPAR2024 <- JASPAR2024::JASPAR2024()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024))
J2024 <- TFBSTools::getMatrixSet(JASPARConnect, opts) #core Jaspar 2024 profiles for plants

J2024_PSBG <- PscanR::ps_build_bg(prom_seq, J2024, BPPARAM = BiocParallel::MulticoreParam(24)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2024_PSBG, "J2024_TAIR9_200u_50d_TAIR.psbg2.txt")