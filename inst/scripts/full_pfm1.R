# Recipe for full_pfms.rds and full_pfm1.rds (not a production background).
# Assisted-by: OpenAI Codex (repair of the historical serialization recipe).
# Rscript full_pfm1.R PROMOTERS.rds OUTPUT_DIRECTORY
# PROMOTERS.rds must contain the first 50 canonical hg38 transcript promoters
# from the original annotation, using -200/+50, before deduplication.
# The bundled objects contain 36 unique sequences and the first 50 CORE
# vertebrate JASPAR2020 motifs. Live annotations can change that count;
# see fixture-provenance.md for the limits of historical reproduction.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !file.exists(args[[2]]))
sequences <- readRDS(args[[1]])
stopifnot(methods::is(sequences, "DNAStringSet"), length(sequences) == 50L)
motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,
    list(collection = "CORE", tax_group = "vertebrates"))[1:50]
full <- PscanR::ps_build_bg(sequences, motifs, fullBG = TRUE,
    BPPARAM = BiocParallel::SerialParam())
dir.create(args[[2]], recursive = TRUE)
saveRDS(full, file.path(args[[2]], "full_pfms.rds"), version = 3)
saveRDS(full[[1]], file.path(args[[2]], "full_pfm1.rds"), version = 3)
