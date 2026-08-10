# PscanR 0.99.0

- Compute the motif enrichment p-value as an exact upper normal tail with `pnorm(z, lower.tail = FALSE)` instead of `1 - pnorm(z)`. The previous formulation saturated at 1.11e-16 and returned exactly 0 for z-scores above roughly 8.3, so the most strongly enriched motifs could not be ranked and their FDR values collapsed to 0. Z-scores, hit scores, positions, strands and oligos are unchanged. This removes the package's only use of `BSDA`, which is no longer a dependency.
- Optimize internal sequence scanning function `.ps_scan_s` using raw byte conversion and vectorized matrix-based scoring.
- Hoist the per-motif forward/reverse-complement score matrices out of `.ps_scan_s` so they are built once per motif and reused across all sequences.
- Add a batched, encode-once scanning fast path: sequences are encoded to integer codes once per analysis (reused across all motifs) and all equal-length sequences are scored against each motif with vectorised matrix operations, with the per-sequence `.ps_scan_s` kept as a fallback for unequal-length input. Results are bitwise identical to the previous implementation.
- Accelerate unambiguous, equal-length sequence sets further with Biostrings' vectorised PWM scoring primitive, while retaining the bounded-memory R scanner for ambiguous or large inputs. Repeated 1000-promoter benchmarks remained exactly identical to the saved reference results.
- Score soft-masked (lower-case) bases as their upper-case equivalent, and return `NA` (instead of `-Inf` or an error) for sequences shorter than the motif or with no scorable window, so degenerate sequences no longer corrupt the foreground/background averages.
- Add regression tests for the sequence scanning kernel.

# PscanR 1.0

- First release of PscanR
