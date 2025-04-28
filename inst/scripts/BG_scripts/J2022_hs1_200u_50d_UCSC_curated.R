gff <- rtracklayer::import('raw_data/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz')
chrominfo <- GenomeInfoDb::getChromInfoFromNCBI("T2T-CHM13v2.0")
GenomeInfoDb::seqlevels(gff) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn)
GenomeInfoDb::seqinfo(gff) <- GenomeInfoDb::Seqinfo(genome="T2T-CHM13v2.0")
txdb <- txdbmaker::makeTxDbFromGRanges(gff, taxonomyId=9606)
GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:24] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
prom_seq <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hs1::BSgenome.Hsapiens.UCSC.hs1, 
                               prom_rng) #promoter sequences

# Remove sequences which identifier starts with XR or XM. Remove also sequences without name

prom_seq <- prom_seq[!grepl('^X[MR]_', names(prom_seq)) & names(prom_seq) != '']

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2022 <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts) #core Jaspar 2022 profiles for vertebrates

J2022_PSBG <- PscanR::ps_build_bg(prom_seq, J2022, BPPARAM = BiocParallel::MulticoreParam(24)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2022_PSBG, "J2022_hs1_200u_50d_UCSC.psbg.txt")
