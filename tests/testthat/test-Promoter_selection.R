test_that("ps_select_promoters ranks select transcripts before fallback", {
  annotation <- data.frame(
    SYMBOL = c(
      "GENE1", "GENE1",
      "GENE2",
      "GENE3", "GENE3",
      "GENE4"
    ),
    REFSEQ = c(
      "NM_000001.1", "NM_000002.1",
      "NM_000003.1",
      "NM_000004.1", "NR_000005.1",
      "XR_000006.1"
    ),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    refseq_id = c("NM_000002.1", "NR_000005.1"),
    refseq_clean = c("NM_000002", "NR_000005"),
    gene_symbol = c("GENE1", "GENE3"),
    selection_source = c("MANE Select", "RefSeq Select"),
    selection_priority = c(1L, 2L),
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    c("GENE1", "GENE2", "GENE3", "GENE4"),
    annotation = annotation,
    select_transcripts = select_tx
  )

  expect_identical(
    selected$refseq_clean,
    c("NM_000002", "NM_000003", "NR_000005", "XR_000006")
  )
  expect_identical(
    selected$selection_source,
    c(
      "MANE Select", "RefSeq NM fallback", "RefSeq Select",
      "RefSeq fallback"
    )
  )
})

test_that("ps_select_promoters can require selected transcripts", {
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE2"),
    REFSEQ = c("NM_000001.1", "NM_000002.1"),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    refseq_clean = "NM_000002",
    selection_source = "MANE Select",
    selection_priority = 1L,
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    c("GENE1", "GENE2"),
    annotation = annotation,
    select_transcripts = select_tx,
    fallback = FALSE
  )

  expect_identical(selected$gene, "GENE2")
  expect_identical(selected$selection_source, "MANE Select")
})

test_that("ps_select_promoters returns selected promoter sequences", {
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE1", "GENE2"),
    REFSEQ = c("NM_000001.1", "NM_000002.1", "NM_000003.1"),
    stringsAsFactors = FALSE
  )
  select_tx <- data.frame(
    refseq_clean = "NM_000002",
    selection_source = "MANE Select",
    selection_priority = 1L,
    stringsAsFactors = FALSE
  )
  promoters <- Biostrings::DNAStringSet(c(
    NM_000001.1 = "AAAA",
    NM_000002.1 = "CCCC",
    NM_000003.1 = "GGGG"
  ))

  selected <- ps_select_promoters(
    c("GENE1", "GENE2"),
    promoter_sequences = promoters,
    annotation = annotation,
    select_transcripts = select_tx,
    return = "sequences"
  )
  mapping <- attr(selected, "ps_gene_promoter_map")

  expect_s4_class(selected, "DNAStringSet")
  expect_identical(names(selected), c("GENE1|NM_000002", "GENE2|NM_000003"))
  expect_identical(unname(as.character(selected)), c("CCCC", "GGGG"))
  expect_identical(mapping$promoter_id, c("NM_000002.1", "NM_000003.1"))
})

test_that("ps_select_promoters can return all mapped transcripts", {
  annotation <- data.frame(
    SYMBOL = c("GENE1", "GENE1", "GENE1"),
    REFSEQ = c("NR_000002.1", "NM_000003.1", "NM_000001.1"),
    stringsAsFactors = FALSE
  )

  selected <- ps_select_promoters(
    "GENE1",
    annotation = annotation,
    mode = "all_transcripts"
  )

  expect_identical(
    selected$refseq_clean,
    c("NM_000001", "NM_000003", "NR_000002")
  )
})
