# PscanR

**PscanR** is an R/Bioconductor package for transcription factor binding motif
enrichment analysis in regulatory DNA sequences. It scans the promoters of a set
of co-regulated or co-expressed genes with a collection of position weight
matrices, usually from JASPAR, and compares the resulting scores against a
background distribution computed over all the promoters of the same organism.
Motifs that score higher in the input set than the background would predict are
reported as candidate common regulators, together with the position of their
predicted binding sites.

It is a more flexible and embeddable implementation of the Pscan algorithm of
Zambelli et al. (2009), which is also available as a web server and a C++
application at <http://www.beaconlab.it/pscan/>. The package builds on
established Bioconductor infrastructure, in particular *Biostrings*,
*TFBSTools* and *BiocParallel*. Its central data structure is the
`PSMatrixList` class, which extends TFBSTools' `PFMatrixList` to carry
background statistics and per-promoter hits alongside the matrices themselves.

## Installation

PscanR is at version 0.99.0 and is **not yet on Bioconductor**. Install it from
GitHub:

```r
pak::pak("Federico77z/PscanR")
```

`pak` resolves Bioconductor repositories on its own, which matters here: PscanR
imports seven Bioconductor packages and suggests fifteen more.

Once the package is accepted into Bioconductor, the usual route will work
instead:

```r
BiocManager::install("PscanR")   # not available yet
```

## Quick start

Everything below ships with the package, so this example runs with no download.
It scans 90 curated human NRF1 target promoters against the JASPAR 2020 core
collection, using the bundled human background for the matching promoter window.

```r
library(PscanR)
library(Biostrings)

# 746 JASPAR 2020 core matrices, and the human background computed for
# promoters 200 bp upstream to 50 bp downstream of the TSS.
J2020 <- readRDS(system.file("extdata", "J2020.rds", package = "PscanR"))
bg_path <- system.file(
    "extdata", "J2020_hg38_200u_50d_UCSC.psbg.txt",
    package = "PscanR"
)
background <- ps_retrieve_bg_from_file(bg_path, J2020)

# The foreground: promoters of transcripts reported as NRF1 targets.
prom_seq <- readRDS(
    system.file("extdata", "prom_seq.rds", package = "PscanR")
)
target <- read.csv(
    system.file("extdata", "nrf1100.txt", package = "PscanR"),
    header = FALSE
)
foreground <- prom_seq[names(prom_seq) %in% target[[1]]]

res <- pscan(foreground, background, BPPARAM = BiocParallel::SerialParam())
head(ps_results_table(res)[, c("NAME", "ZSCORE", "P.VALUE", "FDR")], 5)
```

The scan takes a few seconds on one core and returns:

```
          NAME    ZSCORE      P.VALUE          FDR
MA0506.1  NRF1 17.105626 6.737349e-66 5.026062e-63
MA0615.1 Gmeb1  8.033328 4.743196e-16 1.769212e-13
MA0632.2 TCFL5  7.958568 8.702075e-16 2.163916e-13
MA0641.1  ELF4  6.857577 3.501917e-12 5.733857e-10
MA1483.1  ELF2  6.844280 3.843068e-12 5.733857e-10
```

`ps_results_table()` orders motifs by increasing p-value and reports, for each,
the enrichment z-score, its p-value and the Benjamini-Hochberg FDR. NRF1 comes
first by a wide margin over the second motif, which is the expected answer for a
set of curated NRF1 targets.

## Precomputed backgrounds

A background describes the score distribution of every matrix over all the
promoters of an organism, for one promoter window and one motif collection.
PscanR ships a catalog of precomputed backgrounds covering:

| Species | Assemblies |
| --- | --- |
| *Homo sapiens* | hg38, hs1 |
| *Mus musculus* | mm10, mm39 |
| *Drosophila melanogaster* | dm6 |
| *Arabidopsis thaliana* | TAIR9 |
| *Saccharomyces cerevisiae* | sacCer3 |

each combined with the promoter windows 200 bp upstream to 50 downstream, 450 to
50, 500 to 0, 950 to 50 and 1000 to 0, for the JASPAR 2020, 2022 and 2024 core
collections.

```r
get_availableBG()                # what is available
get_availableBG(details = TRUE)  # the full version-2 catalog

background <- generate_psmatrixlist_from_background(
    "Jaspar2020", "hs", c(-200, 50), "hg38"
)
```

The four arguments are the four things that identify a background: the JASPAR
release, the species, the promoter window relative to the TSS, and the genome
assembly. The assembly is needed only for the two species that have more than
one.

Retrieval goes through ExperimentHub and falls back, with a warning, to the
immutable Zenodo record <https://doi.org/10.5281/zenodo.21821764> when the Hub
cannot be reached. Both routes serve the same archive, and every file is
verified by SHA-256 on arrival. Pass `source = "zenodo"` to go to Zenodo
directly.

**The foreground promoter window must match the background's.** A `950u_50d`
background describes promoters from -950 to +50 relative to the TSS, and a
foreground extracted over any other window is not comparable to it. When no
precomputed background matches the species, window or motif collection an
analysis needs, one can be built with `ps_build_bg()`.

## Documentation

Five vignettes, in the order they are meant to be read:

1. **PscanR Quick Overview**. The exported functions, grouped by task, with a
   one-line description of each.
2. **PscanR Concepts and Input Preparation**. What a background is, how to
   obtain one, and how to get from a list of identifiers to the promoter
   sequences a scan takes.
3. **Motif enrichment analysis with PscanR: mouse retinal promoters**. A
   worked analysis on differentially expressed genes.
4. **Motif enrichment analysis with PscanR: Arabidopsis BPC target promoters**.
   A worked analysis outside the RefSeq identifier scheme.
5. **Building a background for PscanR**. Building, verifying and storing a
   background for an organism or window the catalog does not cover.

```r
browseVignettes("PscanR")
```

## Citation

If PscanR contributes to published work, please cite the paper describing the
algorithm:

> Zambelli F, Pesole G, Pavesi G. Pscan: finding over-represented transcription
> factor binding site motifs in sequences from co-regulated or co-expressed
> genes. *Nucleic Acids Research*, 2009, 37:W247-W252.
> doi:[10.1093/nar/gkp464](https://doi.org/10.1093/nar/gkp464)

`citation("PscanR")` prints the same reference, with a BibTeX entry.

## License

GPL-3.
