# PscanR 0.99.5

- The documentation of `ps_results_table()` and the concepts vignette state
  that `P.VALUE` is one-sided (upper tail, so it tests for enrichment only) and
  that `FDR` adjusts it over all the motifs, before the `FDR` filter is
  applied.
- The documentation of `ps_build_bg_from_table()` states that a
  `PSMatrixList`, such as the output of `ps_retrieve_bg()`, is accepted and
  that its background statistics are replaced by those supplied.

# PscanR 0.99.4

- The README describes installation from Bioconductor, links the
  PscanRBackgrounds package page on Bioconductor, and gives the order in which
  `ps_retrieve_bg()` tries ExperimentHub and the Zenodo records.

# PscanR 0.99.3

- Promoters of different lengths are not comparable, because a promoter's
  best motif score depends on its length. `pscan()`, `pscan_filtered()` and
  `ps_build_bg()` now keep the promoters of the most common length (before,
  those of the longest one), and the warning reports how many promoters were
  kept, their length, and how many were discarded.
- `pscan_filtered()` removes duplicated promoters and promoters of other
  lengths, as `pscan()` does.
- `ps_scan()` warns when it is called directly on promoters of different
  lengths, which it scans one at a time.
- `pscan()`, `pscan_filtered()` and `ps_build_bg()` stop with an error that
  names the motifs when the promoters kept are shorter than any motif in the
  collection. Before, `pscan()` failed with an unrelated z-test error,
  `pscan_filtered()` with "logical subscript contains NAs", and
  `ps_build_bg()` silently stored `NaN` background statistics. `ps_scan()`
  stops with a clear error when no promoter is long enough for its motif.
- `pscan_filtered()` encodes the promoters once for all motifs, as `pscan()`
  does.
- Argument checks are shared by all exported functions and stop with one
  wording, such as "'x' must be a DNAStringSet", without showing internal
  function names. Checks that were missing are added:
  `ps_retrieve_bg_from_file()` requires a file path, `ps_build_bg_from_table()`
  a data.frame,
  `ps_build_bg()` a logical `fullBG`, and `pscan_filtered()` a numeric `n`.
  `ps_write_bg_to_file()` reports an output file or directory it cannot write
  to; the check existed but never ran. `ps_motif_barplot()` accepts an `FDR`
  only between 0 and 1, as `ps_results_table()` does.
- Removed internal code without callers: `.ps_natural_order()`, the stored
  `PSMatrixList` class generator, and the `.PS_PSEUDOCOUNT()` and
  `all_sequences_id()` accessor generics (the latter read a slot that
  `PSMatrix` does not have).
- Tests cover the per-sequence scanner used for promoters of different
  lengths.

# PscanR 0.99.2

- Changes from the Bioconductor review:
  - `ps_motif_barplot()` and the density and hit-score plots use the Okabe-Ito
    colour-blind-safe palette from `grDevices::palette.colors()`.
  - The vignettes locate their bundled files with `system.file()`, so every
    displayed chunk can be run as shown.
  - The funder is listed in `Authors@R`.
  - `PscanRBackgrounds`, the ExperimentHub package with the precomputed
    backgrounds, is in Suggests now that it is part of Bioconductor.
  - Minor code cleanup: no comparisons with `TRUE`, `startsWith()` for prefix
    tests, `toString()` for comma-separated lists.
  - The code follows a lintr configuration (`.lintr`: 4-space indentation,
    80 columns) that CI enforces; the tests were reindented to match.
  - Consistent argument names: `ps_retrieve_bg(jaspar_release = )` (was
    `JASPAR_matrix`) and `ps_scan(use_full_bg = )` (was `use_full_BG`).

# PscanR 0.99.0

- Initial Bioconductor submission.
- `pscan()` implements the Pscan method (Zambelli, Pesole and Pavesi, 2009,
  doi:10.1093/nar/gkp464): it tests a set of promoter sequences for
  transcription factor binding motifs that score higher than in a background
  of all promoters of the organism, and reports z-scores, p-values and FDR.
- Precomputed backgrounds for the JASPAR 2020, 2022 and 2024 CORE collections,
  seven genome assemblies and five promoter windows are retrieved from
  ExperimentHub (package PscanRBackgrounds), or directly from their Zenodo
  records, with `ps_retrieve_bg()` and listed with `ps_available_bg()`.
- Custom backgrounds can be built with `ps_build_bg()` and saved and reloaded
  with `ps_write_bg_to_file()` and `ps_retrieve_bg_from_file()`. Full
  backgrounds keep per-promoter hits, so that `pscan_full_bg()` can test a set
  of transcripts without rescanning sequences.
- `pscan_filtered()` restricts an analysis to motifs that pass a first scan.
- `ps_select_promoters()` and `ps_selection_summary()` choose one promoter per
  gene from a gene list.
- Results are summarised with `ps_results_table()` and visualised with
  `ps_zscore_heatmap()`, `ps_motif_barplot()`, `ps_density_plot()`,
  `ps_hitpos_map()` and `ps_hit_score_plot()`.
- Motif scanning can run in parallel through BiocParallel; it is serial by
  default.
