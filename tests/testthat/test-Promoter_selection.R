test_that("ps_load_select_transcripts returns cached metadata offline", {
  cache_dir <- tempfile("pscanr-select-")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  cache_file <- file.path(
    cache_dir,
    "hg38_mane_refseq_select_transcripts.rds"
  )
  cached_tx <- data.frame(
    transcript_id = c("NM_000546.6", "NM_000492.4"),
    transcript_base = c("NM_000546", "NM_000492"),
    gene_symbol = c("TP53", "CFTR"),
    selection_source = c("MANE Select", "RefSeq Select"),
    selection_priority = c(1L, 2L)
  )
  saveRDS(cached_tx, cache_file)

  local_mocked_bindings(
    .ps_mane_select_transcripts = function() stop("Unexpected network call"),
    .ps_refseq_select_transcripts = function() stop("Unexpected network call"),
    .package = "PscanR"
  )

  expect_identical(
    ps_load_select_transcripts(cache_dir = cache_dir),
    cached_tx
  )
})

test_that("ps_select_promoters ranks select transcripts before fallback", {
  annotation <- data.frame(
    SYMBOL = c(
      "GENE1", "GENE1", "GENE2", "GENE3", "GENE3", "GENE4"
    ),
    REFSEQ = c(
      "NM_000001.1", "NM_000002.1", "NM_000003.1", "NM_000004.1",
      "NR_000005.1", "XR_000006.1"
    ),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    transcript_id = c("NM_000002.1", "NR_000005.1"),
    transcript_base = c("NM_000002", "NR_000005"),
    gene_symbol = c("GENE1", "GENE3"),
    selection_source = c("MANE Select", "RefSeq Select"),
    selection_priority = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    c("GENE1", "GENE2", "GENE3", "GENE4"),
    annotation = annotation,
    select_transcripts = select_tx,
    quiet = TRUE
  )

  expect_identical(
    selected$transcript_base,
    c("NM_000002", "NM_000003", "NR_000005", "XR_000006")
  )
  expect_identical(
    selected$selection_source,
    c("MANE Select", "RefSeq NM", "RefSeq Select", "RefSeq XR")
  )
})

test_that("ps_select_promoters can require selected transcripts", {
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE2"),
    REFSEQ = c("NM_000001.1", "NM_000002.1"),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    transcript_base = "NM_000002",
    selection_source = "MANE Select",
    selection_priority = 1L,
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    c("GENE1", "GENE2"),
    annotation = annotation,
    select_transcripts = select_tx,
    fallback = FALSE,
    quiet = TRUE
  )

  expect_identical(selected$gene, "GENE2")
  expect_identical(selected$selection_source, "MANE Select")
  expect_identical(attr(selected, "ps_unmapped_genes"), "GENE1")
})

test_that("ps_select_promoters returns selected promoter sequences", {
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE1", "GENE2"),
    REFSEQ = c("NM_000001.1", "NM_000002.1", "NM_000003.1"),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    transcript_id = "NM_000002.1",
    transcript_base = "NM_000002",
    gene_symbol = "GENE1",
    selection_source = "MANE Select",
    selection_priority = 1L,
    stringsAsFactors = FALSE
  )
  promoters <- Biostrings::DNAStringSet(c(
    NM_000001.1 = "AAAA", NM_000002.1 = "CCCC", NM_000003.1 = "GGGG"
  ))

  selected <- ps_select_promoters(
    c("GENE1", "GENE2"),
    promoter_sequences = promoters,
    annotation = annotation,
    select_transcripts = select_tx,
    return = "sequences",
    quiet = TRUE
  )

  expect_s4_class(selected, "DNAStringSet")
  # Sequence names now carry the full promoter identifier, not a stripped one.
  expect_identical(names(selected), c("GENE1|NM_000002.1", "GENE2|NM_000003.1"))
  expect_identical(as.character(selected), c("CCCC", "GGGG"), ignore_attr = TRUE)
  expect_identical(
    attr(selected, "ps_gene_promoter_map")$promoter_id,
    c("NM_000002.1", "NM_000003.1")
  )
})

test_that("ps_select_promoters can return all mapped transcripts", {
  annotation <- data.frame(
    SYMBOL = rep("GENE1", 3),
    REFSEQ = c("NR_000002.1", "NM_000003.1", "NM_000001.1"),
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    "GENE1",
    annotation = annotation,
    mode = "all_transcripts",
    quiet = TRUE
  )

  expect_identical(
    selected$transcript_base,
    c("NM_000001", "NM_000003", "NR_000002")
  )
})

# ---------------------------------------------------------------------------
# Transcript identifier schemes
# ---------------------------------------------------------------------------

test_that("transcript schemes strip a version but never a splice variant", {
  expect_identical(
    .ps_transcript_base(c("NM_000546.6", "NR_046018.2"), "refseq"),
    c("NM_000546", "NR_046018")
  )
  # The whole point: TAIR suffixes identify transcripts, so they survive.
  expect_identical(
    .ps_transcript_base(c("AT1G01110.1", "AT1G01110.2"), "tair"),
    c("AT1G01110.1", "AT1G01110.2")
  )
  expect_identical(
    .ps_transcript_base("YAL068W-A", "sgd"), "YAL068W-A"
  )
  # The fallback must not mangle identifiers it does not understand.
  expect_identical(
    .ps_transcript_base("weird.thing.1", "generic"), "weird.thing.1"
  )
})

test_that("scheme detection recognises every supported reference", {
  expect_identical(
    .ps_detect_scheme(c("NM_000546.6", "NR_046018.2", "NP_034436"))$scheme,
    "refseq"
  )
  expect_identical(
    .ps_detect_scheme(c("AT1G01010.1", "AT1G01110.2"))$scheme, "tair"
  )
  expect_identical(
    .ps_detect_scheme(c("YAL069W", "YAL068W-A"))$scheme, "sgd"
  )
  # Mixed grammars are not guessed at; they fall back to the lossless scheme.
  expect_identical(
    .ps_detect_scheme(c("NM_000546.6", "AT1G01010.1"))$scheme, "generic"
  )
  expect_warning(
    .ps_resolve_scheme("auto", promoter_ids = c("NM_000546.6", "AT1G01010.1")),
    "generic"
  )
})

test_that("natural ordering compares embedded numbers numerically", {
  ids <- c("AT4G32850.10", "AT4G32850.2", "AT4G32850.9", "AT4G32850.1")
  expect_identical(
    ids[.ps_natural_order(ids)],
    c("AT4G32850.1", "AT4G32850.2", "AT4G32850.9", "AT4G32850.10")
  )
})

test_that("selection does not depend on the order promoters are supplied", {
  cases <- list(
    refseq = list(
      promoters = c("NM_000001.1", "NM_000002.1", "NR_000003.1",
                    "NM_000004.2"),
      genes = c("G1", "G1", "G1", "G2")
    ),
    tair = list(
      promoters = c("AT1G01110.2", "AT1G01110.1", "AT1G01110.10",
                    "AT1G01160.3"),
      genes = c("AT1G01110", "AT1G01110", "AT1G01110", "AT1G01160")
    ),
    sgd = list(
      promoters = c("YAL069W", "YAL068W-A", "YAL067W-A"),
      genes = c("YAL069W", "YAL068W-A", "YAL067W-A")
    )
  )

  for (nm in names(cases)) {
    case <- cases[[nm]]
    annotation <- data.frame(
      GENE = case$genes, TX = case$promoters, stringsAsFactors = FALSE
    )
    genes <- unique(case$genes)

    forward <- ps_select_promoters(
      genes, promoter_ids = case$promoters, annotation = annotation,
      gene_col = "GENE", transcript_col = "TX", mode = "representative",
      quiet = TRUE
    )
    backwards <- ps_select_promoters(
      genes, promoter_ids = rev(case$promoters),
      annotation = annotation[rev(seq_len(nrow(annotation))), ],
      gene_col = "GENE", transcript_col = "TX", mode = "representative",
      quiet = TRUE
    )

    expect_identical(forward$scheme[[1L]], nm)
    joined <- merge(
      forward[, c("gene", "promoter_id")],
      backwards[, c("gene", "promoter_id")],
      by = "gene"
    )
    expect_identical(nrow(joined), length(genes))
    expect_identical(joined$promoter_id.x, joined$promoter_id.y)
  }
})

test_that("TAIR selection prefers the primary gene model", {
  promoters <- c("AT1G01110.2", "AT1G01110.1", "AT1G01110.10",
                 "AT2G00001.3", "AT2G00001.9")
  annotation <- data.frame(
    GENE = sub("\\..*$", "", promoters), TX = promoters,
    stringsAsFactors = FALSE
  )
  selected <- ps_select_promoters(
    c("AT1G01110", "AT2G00001"),
    promoter_ids = promoters, annotation = annotation,
    gene_col = "GENE", transcript_col = "TX", mode = "representative",
    quiet = TRUE
  )

  # .1 wins when present; otherwise the lowest available variant, with .10
  # ranked after .9 rather than lexicographically before .3.
  expect_identical(
    selected$promoter_id, c("AT1G01110.1", "AT2G00001.3")
  )
  expect_identical(selected$transcript_class, c("primary", "alternative"))
  # No candidate is discarded before ranking.
  expect_identical(selected$n_candidates, c(3L, 2L))
  expect_identical(selected$suffix_role, c("variant", "variant"))
})

test_that("mode = 'select' is refused outside the RefSeq scheme", {
  promoters <- c("AT1G01110.1", "AT1G01110.2")
  annotation <- data.frame(
    GENE = "AT1G01110", TX = "AT1G01110.1", stringsAsFactors = FALSE
  )
  expect_error(
    ps_select_promoters(
      "AT1G01110", promoter_ids = promoters, annotation = annotation,
      gene_col = "GENE", transcript_col = "TX", mode = "select", quiet = TRUE
    ),
    "MANE Select"
  )
})

test_that("protein accessions are excluded and reported", {
  promoters <- c("NM_000001.1", "NM_000002.1")
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE1", "GENE2"),
    REFSEQ = c("NM_000001", "NP_000001", "NM_000002"),
    stringsAsFactors = FALSE
  )
  selected <- ps_select_promoters(
    c("GENE1", "GENE2"), promoter_ids = promoters, annotation = annotation,
    mode = "representative", quiet = TRUE
  )
  expect_false(any(grepl("^NP_", selected$promoter_id)))
  summary <- ps_selection_summary(selected)
  expect_identical(summary$overall$candidates_excluded_protein, 1L)
  expect_identical(summary$overall$genes_mapped, 2L)
})

test_that("ps_selection_summary reports unmapped genes and tie-breaks", {
  promoters <- c("NM_000001.1", "NM_000002.1")
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE1"),
    REFSEQ = c("NM_000001", "NM_000002"),
    stringsAsFactors = FALSE
  )
  selected <- ps_select_promoters(
    c("GENE1", "GENE_ABSENT"), promoter_ids = promoters,
    annotation = annotation, mode = "representative", quiet = TRUE
  )
  summary <- ps_selection_summary(selected)
  expect_identical(summary$overall$genes_requested, 2L)
  expect_identical(summary$overall$genes_unmapped, 1L)
  expect_identical(attr(selected, "ps_unmapped_genes"), "GENE_ABSENT")
  # Two equally-ranked NM_ candidates: honest about how it was decided.
  expect_identical(selected$decided_by, "tie_break")
  expect_identical(summary$overall$genes_decided_by_tie_break, 1L)
})
