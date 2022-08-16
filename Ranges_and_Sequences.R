library(GenomicFeatures)
library(GenomicRanges)

txdb <- makeTxDbFromUCSC(genome="sacCer2", tablename="ensGene")

txdb <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene",
                 transcript_ids=NULL,
                 circ_seqs=NULL,
                 url="http://genome.ucsc.edu/cgi-bin/",
                 goldenPath.url="https://hgdownload.soe.ucsc.edu/",
                 taxonomyId=NA,
                 miRBaseBuild=NA)