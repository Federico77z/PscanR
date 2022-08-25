library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)


txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") #import gtf annotation from ucsc
seqlevels(txdb) <- seqlevels(txdb)[1:24] #use only annotation on canonical chromosomes

export(txdb, "./Annotation/hg38.ucsc.bed", format = "BED")  #export
export(txdb, "./Annotation/hg38.ucsc.bed", format = "GTF")

#trnames <- transcripts(txdb, columns = "TXNAME")$TXNAME
#trasbed <- asBED(txdb)

#txstrand <- select(txdb, trnames, columns = c("TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND", "TXTYPE"), keytype = "TXNAME")$TXSTRAND
#txstrand_minus <- which(txstrand == "-")
#txstrand_plus <- which(txstrand == "+")

#select(txdb, trnames[txstrand_minus], columns = c("TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND", "TXTYPE"), keytype = "TXNAME", )

prom_rng <- promoters(txdb, upstream = 500, downstream = 0, use.names = TRUE) #Define promoters

#prom_seq <- getPromoterSeq(asBED(txdb), BSgenome.Hsapiens.UCSC.hg38, upstream = 500, downstream = 0) #get promoter sequences
#prom_seq <- extractUpstreamSeqs(BSgenome.Hsapiens.UCSC.hg38, txdb, width = 500)
#names(prom_seq) <- keys(txdb, "TXNAME")

prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, prom_rng)  #Get promoter sequences
prom_seq_unique <- unique(prom_seq)

#tx_plus_pos <- which(as.character(strand(prom_rng)) == "+")
#tx_minus_pos <- which(as.character(strand(prom_rng)) == "-")

writeXStringSet(prom_seq, "./FASTA/hg38_prom_500_0.fasta") #Write multifasta file
writeXStringSet(prom_seq_unique, "./FASTA/hg38_prom_500_0.unique.fasta")

pcharSpl <- strsplit(as.character(prom_seq_unique$NR_046018.2), "")[[1]]

pcharSpl10 <- pcharSpl[1:10]

p_alpha <- c(1L,2L,3L,4L)
names(p_alpha) <- c("A", "C", "G", "T")
pwm_idx <- cbind(p_alpha[pcharSpl10], 1:10)
pwm[pwm_idx]

substring(pchar, 1:(nchar(pchar) - 10 + 1), 10:nchar(pchar))

