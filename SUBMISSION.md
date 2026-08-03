# Submitting PscanR (and PscanRBackgrounds) to Bioconductor

PscanR is submitted together with its companion **ExperimentHub** data package
**PscanRBackgrounds**, which distributes the precomputed motif background models
that used to be downloaded from a personal GitHub repository. The two are
reviewed together as a *package stack*.

This file documents the steps a maintainer performs by hand (account setup,
data upload, the Contributions issue, review). The code changes themselves are
already in the two packages.

---

## 0. One-time maintainer prerequisites

- Create an account on the Bioconductor support site
  <https://support.bioconductor.org> using **the same e-mail as the `Maintainer`
  field** (`federico.zambelli@unimi.it`).
- Subscribe to the **bioc-devel** mailing list
  <https://stat.ethz.ch/mailman/listinfo/bioc-devel>.
- Develop/test against **Bioconductor devel** (currently 3.23) on the matching
  R version. Install the check tools once:
  ```r
  BiocManager::install(c("BiocCheck", "ExperimentHub", "ExperimentHubData",
                         "GenomeInfoDbData"))
  ```

---

## 1. PscanRBackgrounds — register and upload the background data

The version-2 background tables are **not** shipped inside PscanR. They are
distributed as a single ZIP archive, `PscanR_backgrounds_v2.zip`, hosted on the
Bioconductor ExperimentHub S3 bucket and described by
`inst/extdata/metadata.csv` in the data package. PscanR locates it by
`preparerclass == "PscanRBackgrounds"` and `title == "PscanR_backgrounds_v2"`,
so no `EH*****` id is hard-coded.

The archive contents are indexed by the catalog PscanR bundles at
`inst/extdata/PscanR_background_catalog_v2.tsv` (105 validated combinations).
Keep that file synchronised with `catalog.tsv` in the background repository.

1. Regenerate / validate the metadata table:
   ```bash
   cd PscanRBackgrounds
   Rscript inst/scripts/make-metadata.R            # writes inst/extdata/metadata.csv
   ```
   ```r
   ExperimentHubData::makeExperimentHubMetadata("PscanRBackgrounds")  # must pass cleanly
   ```
2. **Upload the data to S3.** Ask the Bioconductor Hubs team for write access to
   the staging bucket by e-mailing **hubs@bioconductor.org** (mention the package
   name `PscanRBackgrounds`). Then upload `PscanR_backgrounds_v2.zip` under the
   prefix `PscanRBackgrounds/` (matching the `RDataPath` column), following the
   current "Upload to S3" instructions in the *Creating ExperimentHub packages*
   vignette (`vignette("CreateAnExperimentHubPackage", package = "ExperimentHubData")`).
3. Push `PscanRBackgrounds` to a public GitHub repo with `Version: 0.99.0`.

> Note: the generation and publishing pipeline for the archive lives in the
> `PscanR_backgrounds` repository (`scripts/backgrounds.R`), which also produces
> the SHA-256 recorded in both the catalog and `R/Helper_functions.R`.

---

## 2. Local checks — both packages must be clean

On Bioconductor devel, for each package:

```bash
R CMD build PscanRBackgrounds
R CMD BiocCheck PscanRBackgrounds_0.99.0.tar.gz

R CMD build PscanR
R CMD check  PscanR_0.99.0.tar.gz
R CMD BiocCheck PscanR_0.99.0.tar.gz
```

Target: **0 ERROR / 0 WARNING** from both `R CMD check` and `BiocCheck`.

PscanR's examples and tests run offline against the bundled example background
in `inst/extdata/`, so `R CMD check` does not need Hub access. Building the
vignettes that call `generate_psmatrixlist_from_background()` does require the
archive to be reachable; before the data is uploaded (step 1), retrieval falls
back to the immutable Zenodo record (<https://doi.org/10.5281/zenodo.21821764>),
so the vignettes still build.

---

## 3. Open the submission issue

1. Go to <https://github.com/Bioconductor/Contributions/issues> and open a **New
   issue** with the PscanR repository URL as described in the Contributions
   `README`.
2. In the issue body, add the companion package so both are built together:
   ```
   AdditionalPackage: https://github.com/Federico77z/PscanRBackgrounds
   ```
3. Confirm you have read the contribution guidelines; the webhook adds the
   `Single Package Builder` (SPB) which builds both packages on Linux, macOS and
   Windows on every push.

---

## 4. Review loop

- Fix anything the SPB reports until the build is green on all three platforms.
- A Bioconductor reviewer is assigned; the **Hubs team** finalises the S3 data
  and assigns the `EH*****` id for PscanRBackgrounds.
- Respond to review comments by pushing new commits (bump the `z` in `0.99.z`
  for substantive changes).

---

## 5. After acceptance

- Both packages are added to Bioconductor **devel** and released at the next
  Bioconductor release cycle.
- At the first release the versions move to `1.0.0` automatically.
- From then on, the backgrounds are updated by publishing a new archive version
  from the `PscanR_backgrounds` repository and editing `PscanRBackgrounds`
  (`metadata.csv` + a new `SourceVersion`/`Title`) — independently of the PscanR
  software release, which is exactly why ExperimentHub was chosen over bundling
  the data.

---

## How PscanR retrieves backgrounds

- `generate_psmatrixlist_from_background()` and `get_availableBG()` take a
  `source` argument defaulting to `"experimenthub"`.
- `"experimenthub"` resolves the archive through `ExperimentHub` and falls back
  automatically to the same immutable Zenodo record when the Hub is unavailable;
  `"zenodo"` targets that record directly. Both verify the archive SHA-256 and
  cache it with `BiocFileCache`.
- `"github"` is retained as an explicit legacy backend and is the only source
  that exposes version-1 resources.
- The `httr` dependency and its insecure `httr::set_config(ssl_verifypeer = 0L)`
  override were removed; the legacy path now uses `utils::download.file()`.
