source("load_packages.R", local = TRUE, echo = FALSE)

source("R/PSMatrix_class.R", local = TRUE, echo = FALSE)
source("R/Build_Background.R", local = TRUE, echo = FALSE)
source("R/Helper_functions.R", local = TRUE, echo = FALSE)
source("R/Scan_and_post_processing.R", local = TRUE, echo = FALSE)


#txdb <- makeTxDbFromUCSC(genome="mm10", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene") #import gtf annotation from UCSC
seqlevels(txdb) <- seqlevels(txdb)[1:24] #use only annotations on canonical chromosomes

prom_rng <- promoters(txdb, upstream = 450, downstream = 50, use.names = TRUE) 
prom_seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm10, prom_rng) #promoter sequences

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2020 <- getMatrixSet(JASPAR2020, opts) #core Jaspar 2020 profiles for vertebrates

J2020_PSBG <- ps_build_bg(prom_seq, J2020, BPPARAM = MulticoreParam(12)) #Build Pscan Background

ps_write_bg_to_file(J2020_PSBG, "J2020_mm10_450u_50d_UCSC_refGene.psbg.txt")

rm(list = c("J2020", "opts", "prom_seq", "prom_rng"))
