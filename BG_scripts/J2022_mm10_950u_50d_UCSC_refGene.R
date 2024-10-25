source("load_packages.R", local = TRUE, echo = FALSE)

source("../R/PSMatrix_class.R", local = TRUE, echo = FALSE)
source("../R/Build_Background.R", local = TRUE, echo = FALSE)
source("../R/Helper_functions.R", local = TRUE, echo = FALSE)
source("../R/Scan_and_post_processing.R", local = TRUE, echo = FALSE)

#txdb <- makeTxDbFromUCSC(genome="mm10", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene") #import gtf annotation from UCSC
seqlevels(txdb) <- seqlevels(txdb)[1:24] #use only annotations on canonical chromosomes

prom_rng <- promoters(txdb, upstream = 950, downstream = 50, use.names = TRUE) 
prom_seq <- getSeq(x = BSgenome.Mmusculus.UCSC.mm10, prom_rng) #promoter sequences
export(prom_rng, con = "mm10_950_50.gtf", format = "gtf")

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2022 <- getMatrixSet(JASPAR2022, opts) #core Jaspar 2022 profiles for vertebrates

J2022_PSBG <- ps_build_bg(prom_seq, J2022, BPPARAM = MulticoreParam(36)) #Build Pscan Background

ps_write_bg_to_file(J2022_PSBG, "J2022_mm10_950u_50d_UCSC_refGene.psbg.txt")

rm(list = c("J2022", "opts", "prom_seq", "prom_rng"))
