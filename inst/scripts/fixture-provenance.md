# Example-data provenance and reproduction

The files below are examples, not newly estimated production backgrounds.
`../DATA_SOURCES.md` records attribution and upstream reuse terms. The package
contains serialized R objects (`readRDS`), text tables (`read.delim`), and a
versioned transcript list (`readLines`). Do not overwrite installed fixtures
when reproducing them. Use a new temporary output directory and retain
`sessionInfo()`, input checksums and the PscanR commit with the results.

## Human introductory fixtures

- `J2020.rds`: 746 JASPAR 2020 CORE vertebrate PFMs in a TFBSTools PFMatrixList.
  `J2020.R` reconstructs this object using `getMatrixSet()` and `saveRDS()`.
- `J2020_hg38_200u_50d_UCSC.psbg.txt`: historical short background for the same
  motifs, using unique canonical hg38 ncbiRefSeqCurated transcript promoters,
  -200/+50 relative to each strand-aware TSS. The accompanying
  `J2020_hg38_200u_50d_UCSC_curated.R` describes extraction and scanning.
  This historical teaching background is not the v2 catalog entry.
- `nrf1100.txt`: historical human NRF1 target-transcript example from the
  original Pscan web interface, retained in PscanR with the RefSeq versions
  used by its promoter fixture. This origin was confirmed by the maintainer
  on 2026-09-15. Cite Zambelli, Pesole and Pavesi (2009),
  https://doi.org/10.1093/nar/gkp464, for Pscan and its web interface.
  The original experiment/publication underlying this particular list is
  unknown. The Pscan paper is a software citation, not evidence identifying
  that experiment. Use the list as a historical software demonstration.
- `prom_seq.rds`: DNAStringSet of the 90 matching NRF1 target promoters.
  `prom_seq.R` extracts -200/+50 promoters and matches complete RefSeq names.
  The original UCSC snapshot is not preserved for these historical fixtures;
  a current UCSC query may not reproduce their coordinates. Do not substitute
  versionless matches and call the result the same fixture.
- `pfm1.rds`: MA0506.1 scanned against those 90 promoters with the bundled
  background. `pfm1.R` recomputes this analysis entirely from installed inputs.
  The retained historical object differs from the current scanner in 14 hit
  positions/strands/oligos; score and enrichment slots agree under `all.equal`.
  It is a historical accessor/plotting example, not an exact scanner reference.
  The original scanner revision was not recorded. Compare recomputed output
  explicitly before choosing to replace this ordinary teaching fixture.
- `full_pfms.rds`, `full_pfm1.rds`: first 50 motifs / first motif from a full
  background on the first 50 canonical hg38 promoters selected by the original
  recipe (36 unique sequences). `full_pfm1.R` now uses `saveRDS()` for both
  files and the correct plural list filename. Supply the original promoter
  inputs to reproduce historical values; they are not recoverable from the
  retained hit oligos alone.

The original annotation-dependent recipes document the operations, not an
unverifiable promise of bit-for-bit recreation from changing web services.

## Catalog

`PscanR_background_catalog_v2.tsv` is the subset with `background_version == 2`
of the generation repository's `catalog.tsv`, with `BG_files/` in the artifact
column replaced by `backgrounds/`, the path inside the immutable ZIP.
The 105 rows, all other columns, and checksums must agree. The source release is
https://doi.org/10.5281/zenodo.21821764 and its preparation code is at
https://github.com/Federico77z/PscanRBackgrounds/tree/7516eee.

## Mouse retinal analysis (`vignettes/mouse_hd_retina/`)

Source: Cano-Cano et al. (2024), Scientific Reports 14:4176,
https://doi.org/10.1038/s41598-024-54347-8, Supplementary Table S2 in
`41598_2024_54347_MOESM1_ESM.pdf`.

1. Extract the two DESeq2 tables (retina and striatum) from the PDF, preserving
   Ensembl gene IDs and fold-change signs. The local preparation used
   Ghostscript `-sDEVICE=txtwrite`. Validate retina down/up counts 1078/575 and
   striatum down/up counts 763/176. Split by direction and tissue intersection:
   a = retina-only down, b = shared down, c = striatum-only down,
   d = retina-only up, e = shared up, f = striatum-only up.
   Write long-form IDs/subset membership to `paper_subset_membership.tsv`.
2. Load the generation repository's v2 mm10 -950/+50 annotation snapshot;
   `mm10_950u_50d_refseq_snapshot.rds` is its bundled copy. Extract promoters
   with BSgenome.Mmusculus.UCSC.mm10. Map Ensembl genes to RefSeq using
   org.Mm.eg.db; preserve all candidate mappings before selection.
3. Run `ps_select_promoters(mode="representative")` with those sequences and
   the ENSEMBL/REFSEQ annotation columns, then `pscan()` on each subset against
   JASPAR2020/mm10/-950/+50 background v2. The original full-collection runs
   used 746 motifs. The published background window agrees with this example.
4. `promoter_selection_example.rds` stores ten genes, candidate annotation
   rows and their real promoter DNAStringSet. It illustrates RefSeq naming
   and selection; it is not a random or genome-wide background.
5. `uc2_vignette_data.rds` stores the six-set QC, timings and IRF/STAT summaries,
   complete enrichment tables for sets d/f, the d-set promoter mapping, and
   PSMatrixList subsets with motifs MA1418.1, MA0652.1, MA0517.1, MA1623.1,
   MA1513.1, MA1650.1, MA0162.4 and MA0668.1. Its filtered_irf3 component is
   `pscan_filtered()` on set d, anchored on MA1418.1, n=1, using all 746 motifs.
   Recompute enrichment before selecting the plotting motifs, so saved FDR
   values still refer to the complete tested collection. Serialize the list
   with `saveRDS()`. Timings are historical measurements, not reproducible
   numerical constants.

The complete original scripts are currently maintained in the local UseCase
workspace. These installed recipes describe the necessary transformations;
a public immutable archive of the exact full-analysis scripts and input
mapping versions would strengthen reproducibility for review.

## Arabidopsis BPC analysis (`vignettes/arabidopsis_bpc/`)

Source: Caselli et al. (2025), Plant Molecular Biology 116:4,
https://doi.org/10.1007/s11103-025-01662-x. Supplementary Table 1 supplies
1247 downregulated and 1006 upregulated genes; Tables 3 and 4 supply published
motif results and DAP-seq overlaps.

1. Extract TAIR gene IDs and expression direction from Table 1 to the two
   columns of `arabidopsis_bpc_gene_sets.tsv`. Preserve the source ordering.
2. Extract nuclear TAIR9 promoters from the background release's annotation
   snapshots and BSgenome.Athaliana.TAIR.TAIR9. Use the 1000/0 and 950/50
   windows separately, with their matching JASPAR2020 plant backgrounds.
3. Map TAIR transcript IDs to gene IDs by removing the final splice suffix in
   the annotation table only. Supply the complete transcript names to
   `ps_select_promoters(scheme="tair", mode="representative")`; compare with
   the documented primary/all-transcript analyses. Scan each foreground
   against all 530 plant motifs before selecting plotting subsets.
4. `arabidopsis_bpc_results.rds` contains full enrichment tables, mapping and
   coverage summaries, published-result comparisons, motif-family summaries,
   timings and selected PSMatrixList hits. Its main analysis is the 1000/0
   representative-promoter downregulated set. `filtered_bpc6` uses
   pscan_filtered with the MA1402.1 anchor, n=1; `filtered_subset` preserves
   MA1402.1, MA1381.1, MA1372.1, MA1274.1, MA0931.1 and MA1081.1 hits.
   Keep full-collection adjusted p-values when subsetting the plotted motifs.
5. `promoter_selection_example.rds` is generated by
   `make-arabidopsis-example.R`: the first twelve sorted downregulated genes
   with annotated promoters, all their candidate transcript sequences,
   coordinates and source hashes. It is a new, fully reproducible small
   example, independent of the saved whole-analysis statistics.

The paper used TAIR10 and -1000/+100, so the TAIR9 examples are an adaptation,
not an exact reconstruction of the original promoter universe. Archive the
full-analysis preparation scripts and their environment before submission.

## Custom yeast backgrounds (`vignettes/custom_background/`)

1. Extract sacCer3 nuclear transcript promoters from the versioned sgdGene
   snapshot using -950/+50 and BSgenome.Scerevisiae.UCSC.sacCer3. Retain all
   original transcript names and duplicate sequences. Apply bounds/width/N
   exclusions consistently with the generation pipeline.
2. Obtain JASPAR2020 CORE fungi motifs and run `ps_build_bg()`. Write the short
   table to `J2020_sacCer3_950u_50d_UCSC.rebuild.txt`. Compare its motif moments
   with the archived v2 entry. The vignette records the historical comparison.
3. `fullbg_motifs.tsv` lists the fifteen IDs, names and reasons for their
   selection. Build a full background using exactly those motifs and all
   original promoter names; `saveRDS(..., compress="xz")` produces
   `sacCer3_950u_50d_fullBG_15motifs.rds`.
4. Select standard yeast gene names matching `^RP[LS][0-9]` from
   org.Sc.sgd.db's GENENAME map, map to SGD systematic IDs, sort by gene name,
   and retain IDs in the background. Rebuild the DNAStringSet from characters
   before serializing `rp_promoters.rds` with xz compression, avoiding retention
   of the full shared sequence pool.
5. `timings_1core.tsv` and `timings_8core.tsv` are direct wall-clock measurements
   of separate serial/eight-worker runs. `ucbg_summary.rds` stores validation,
   sizes, input counts and the foreground comparison. To reproduce substantive
   results, compare `pscan(rp_promoters, background)` with
   `pscan_fullBG(names(rp_promoters), full_background)`. Hardware-dependent
   timings and compressed sizes are not acceptance criteria for equality.

The new vignette's small round trip is run directly from these real input
sequences; it never derives benchmark timings or substitutes for production
backgrounds. 
