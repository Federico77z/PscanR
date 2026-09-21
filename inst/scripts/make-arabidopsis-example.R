# Assisted-by: OpenAI Codex.
# Reproduce the small TAIR9 promoter-selection example from the pinned
# annotation snapshot in the generation repository at commit 7516eee.
# Rscript make-arabidopsis-example.R /path/to/PscanRBackgrounds OUTPUT.rds
# The output must not exist. Requires PscanR and BSgenome.Athaliana.TAIR.TAIR9.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !file.exists(args[[2]]))
repo <- normalizePath(args[[1]], mustWork = TRUE)
catalog <- utils::read.delim(file.path(repo, "catalog.tsv"))
entry <- subset(catalog, assembly == "TAIR9" & upstream == 1000 &
    downstream == 0 & jaspar_release == 2020 & background_version == 2)
stopifnot(nrow(entry) == 1L)
snapshot <- readRDS(file.path(repo, entry$annotation_snapshot))
gene_sets <- utils::read.delim(system.file("extdata", "vignettes",
    "arabidopsis_bpc", "arabidopsis_bpc_gene_sets.tsv", package = "PscanR"))
tx <- snapshot$transcripts
# See fixture-provenance.md for the publication and its supplementary table.
genes <- head(sort(intersect(gene_sets$gene_id[gene_sets$direction == "down"],
    sub("\\.[0-9]+$", "", tx$transcript_id))), 12L)
tx <- tx[sub("\\.[0-9]+$", "", tx$transcript_id) %in% genes, ]
anchor <- GenomicRanges::GRanges(tx$seqname,
    IRanges::IRanges(tx$tss, width = 1L), strand = tx$strand)
ranges <- GenomicRanges::promoters(anchor, upstream = 1000, downstream = 0)
genome <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9
sequences <- Biostrings::getSeq(genome, ranges)
names(sequences) <- tx$transcript_id
example <- list(genes = genes,
    annotation = data.frame(GENE = tx$gene_id, TX = tx$transcript_id),
    promoter_sequences = sequences, coordinates = tx,
    source = entry[c("assembly", "upstream", "downstream", "annotation_hash",
        "promoter_hash", "annotation_snapshot")])
saveRDS(example, args[[2]], compress = "xz", version = 3)

