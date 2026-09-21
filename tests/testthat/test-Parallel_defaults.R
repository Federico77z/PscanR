# Assisted-by: OpenAI Codex. Defaults must remain explicitly serial.
test_that("default analysis stays serial and agrees with explicit serial runs", {
    for (fun in list(pscan, ps_build_bg, pscan_filtered)) {
        param <- eval(formals(fun)$BPPARAM)
        expect_s4_class(param, "SerialParam")
        expect_equal(BiocParallel::bpworkers(param), 1L)
    }
    motifs <- readRDS(system.file("extdata", "J2020.rds", package="PscanR"))[1:2]
    sequences <- readRDS(system.file("extdata", "prom_seq.rds", package="PscanR"))
    full <- ps_build_bg(sequences, motifs, fullBG=TRUE)
    expect_identical(full, ps_build_bg(sequences, motifs, fullBG=TRUE,
        BPPARAM=BiocParallel::SerialParam()))
    foreground <- sequences[1:5]
    expect_identical(pscan(foreground, full), pscan(foreground, full,
        BPPARAM=BiocParallel::SerialParam()))
    # Stored-hit retrieval does not scan or accept a parallel backend.
    expect_identical(ps_hits_table(pscan_fullBG(names(foreground), full)[[1]]),
        ps_hits_table(pscan(foreground, full)[[1]]))
})
