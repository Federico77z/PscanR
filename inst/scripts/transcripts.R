gtf <- rtracklayer::import("/home/dbetelli/PscanR/test_and_scripts/prove_ENSEMBL/gencode.v47.basic.annotation.gff3.gz")

genes <- c('ENSG00000171723.16', 'ENSG00000112715.26', 'ENSG00000142156.16', 'ENSG00000170873.19')

transcripts <- gtf[gtf$gene_id %in% genes & !is.na(gtf$transcript_support_level) & gtf$transcript_support_level == '1', ]$transcript_id

transcripts <- unique(transcripts)

saveRDS(transcripts, '/home/dbetelli/PscanR/inst/extdata/transcripts.rds')
