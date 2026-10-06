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
