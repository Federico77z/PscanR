library("GenomicFeatures")
library("GenomicRanges")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("TFBSTools")
library("JASPAR2022")

source("../PSMatrix_class.R", local = TRUE, echo = FALSE)
source("../Build_Background.R", local = TRUE, echo = FALSE)
source("../Helper_functions.R", local = TRUE, echo = FALSE)
source("../Scan_and_post_processing.R", local = TRUE, echo = FALSE)


txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
seqlevels(txdb) <- seqlevels(txdb)[1:24] #use only annotations on canonical chromosomes

prom_rng <- promoters(txdb, upstream = 450, downstream = 50, use.names = TRUE) 
prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, prom_rng) #promoter sequences

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates"

J2022 <- getMatrixSet(JASPAR2022, opts) #core Jaspar 2022 profiles for vertebrates

J2022_PSBG <- ps_build_bg(prom_seq, J2022, BPPARAM = MulticoreParam(36)) #Build Pscan Background

ps_write_bg_to_file(J2022_PSBG, "J2022_hg38_450u_50d_UCSC.psbg.txt")

rm(list = c("J2022", "opts", "prom_seq", "prom_rng"))
