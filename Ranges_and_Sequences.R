library(GenomicFeatures)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)


txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") #import gtf annotation from ucsc
seqlevels(txdb) <- seqlevels(txdb)[1:24] #use only annotation on canonical chromosomes

export(txdb, "./Annotation/hg38.ucsc.bed", format = "BED")  #export
export(txdb, "./Annotation/hg38.ucsc.bed", format = "GTF")

prom_seq <- getPromoterSeq(asBED(txdb), BSgenome.Hsapiens.UCSC.hg38, upstream = 500, downstream = 0) #get promoter sequences
#prom_seq <- extractUpstreamSeqs(BSgenome.Hsapiens.UCSC.hg38, txdb, width = 500)
names(prom_seq) <- keys(txdb, "TXNAME")
writeXStringSet(prom_seq, "./FASTA/hg38_prom_500_0.fasta")
