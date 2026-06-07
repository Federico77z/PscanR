# PscanR 0.99.0

- Optimize internal sequence scanning function `.ps_scan_s` using raw byte conversion and vectorized matrix-based scoring.
- Hoist the per-motif forward/reverse-complement score matrices out of `.ps_scan_s` so they are built once per motif and reused across all sequences.
- Add a batched, encode-once scanning fast path: sequences are encoded to integer codes once per analysis (reused across all motifs) and all equal-length sequences are scored against each motif with vectorised matrix operations, with the per-sequence `.ps_scan_s` kept as a fallback for unequal-length input. Results are bitwise identical to the previous implementation.
- Score soft-masked (lower-case) bases as their upper-case equivalent, and return `NA` (instead of `-Inf` or an error) for sequences shorter than the motif or with no scorable window, so degenerate sequences no longer corrupt the foreground/background averages.
- Add regression tests for the sequence scanning kernel.

# PscanR 1.0

- First release of PscanR