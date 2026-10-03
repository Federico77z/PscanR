# Data sources and attribution

PscanR code is distributed under the license in DESCRIPTION. Upstream data
retain their own attribution and reuse conditions; the package license does
not replace them. Preparation and object-by-object descriptions are installed
in `scripts/fixture-provenance.md`.

| Bundled material | Source and attribution | Reuse information |
| --- | --- | --- |
| J2020.rds, motifs in PSMatrix objects and motif-derived score tables | JASPAR CORE 2020; Fornes et al., doi:10.1093/nar/gkz1001 | JASPAR website identifies CC BY 4.0; cite the database publication and retain motif identifiers. See https://jaspar.elixir.no/faq/. |
| Version-2 background catalog | PscanR background release 2, https://doi.org/10.5281/zenodo.21821764 | CC BY 4.0; credit Federico Zambelli and Giulio Pavesi. |
| Mouse retinal gene lists and summaries | Cano-Cano et al. (2024), Scientific Reports 14:4176, doi:10.1038/s41598-024-54347-8, supplementary Table S2 | Article CC BY 4.0, subject to exceptions in individual material credit lines. Tables are extracted, regrouped and analyzed; these transformations are described in the recipe. |
| Arabidopsis gene lists and summaries | Caselli et al., Plant Molecular Biology 116:4 (2026; online December 2025), doi:10.1007/s11103-025-01662-x, supplementary Tables 1, 3 and 4 | Article CC BY 4.0, subject to exceptions in individual material credit lines. The analysis adapts the promoter definition and annotation, as documented in the vignette. |
| Promoter sequences and transcript mappings | UCSC hg38/mm10/sacCer3, NCBI RefSeq, SGD and TAIR9; corresponding BSgenome and annotation packages | Cite the assembly and annotation providers. Upstream source terms must be reviewed for the particular archived inputs; especially confirm redistribution terms for historical TAIR9 inputs. No additional upstream license is asserted here. |
| nrf1100.txt and derived human teaching examples | Historical NRF1 example from the original Pscan web interface, as confirmed by the maintainer; Zambelli, Pesole and Pavesi (2009), doi:10.1093/nar/gkp464 | The paper documents the software/interface. The original experiment underlying this list is unknown; these fixtures demonstrate software operation. Historical annotation-snapshot limitations are described in the preparation recipe. |

CC BY 4.0: https://creativecommons.org/licenses/by/4.0/.
JASPAR 2022: doi:10.1093/nar/gkab1113; JASPAR 2024: doi:10.1093/nar/gkad1059.
The immutable production archive contains computed statistics, not raw genomes.
