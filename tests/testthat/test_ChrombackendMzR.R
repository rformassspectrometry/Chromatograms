test_that("ChromBackendMzR works", {
    expect_true(isReadOnly(be_mzr))
    expect_false(be_mzr@inMemory)
    expect_false(identical(peaksData(be_mzr),.peaksData(be_mzr)))
    expect_true(identical(length(peaksData(be_mzr)), length(.peaksData(be_mzr))))
    expect_true(identical(
        chromData(be_mzr),
        fillCoreChromVariables(.chromData(be_mzr))
    ))
    expect_false(supportsSetBackend(be_mzr))
})

test_that("backendInitialize works", {
    expect_equal(
        backendInitialize(ChromBackendMzR(), files = character()),
        ChromBackendMzR()
    )
    expect_error(
        backendInitialize(ChromBackendMzR(), files = numeric(2)),
        "must be"
    )
})

test_that("show method for ChromBackendMzR works correctly", {
    expect_output(show(be_mzr), "file")
    expect_output(show(be_mzr), "138 chromatograms")
    expect_output(show(be_mzr), "chromIndex")
    expect_output(show(be_mzr), "msLevel")
    expect_output(show(be_mzr), "mz")
    tmp <- be_mzr
    tmp@inMemory <- TRUE
    expect_output(show(tmp), "Peaks data is cached in memory")
})

test_that("replacement method works", {
    tmp <- be_mzr
    peaksData(tmp)[[1]] <- peaksData(be_mzr)[[1]] + 1
    expect_false(identical(peaksData(tmp)[[1]], peaksData(be_mzr)[[1]]))
    expect_true(identical(peaksData(tmp)[[1]], peaksData(be_mzr)[[1]] + 1))
    expect_true(.inMemory(tmp))
    expect_false(.inMemory(be_mzr))

    cd <- chromData(tmp)
    cd$mz <- 1
    expect_true(!identical(chromData(tmp), cd))
    chromData(tmp) <- cd
    expect_equal(chromData(tmp), cd)
})

test_that("backendParallelFactor works", {
    expect_equal(levels(backendParallelFactor(be_mzr)), unique(dataOrigin(be_mzr)))
})

test_that("error message work", {
    expect_error(
        peaksData(be_mzr, columns = "notacolumn"),
        "requested peaks variables"
    )
})

test_that("chromExtract works correctly for ChromBackendMzR", {
    cd <- chromData(be_mzr)
    first_cd <- cd[100:102,]

    ## We'll use known retention time ranges and identifiers from chromData
    peak_tbl <- data.frame(
        rtMin = c(13.7, 16.2, 18 ),
        rtMax = c(18.1,19,20.6 ),
        dataOrigin = first_cd$dataOrigin,
        chromIndex = c(100, 101, 102)
    )
    out <- chromExtract(be_mzr, peak_tbl, by = c("dataOrigin", "chromIndex"))

    expect_s4_class(out, "ChromBackendMzR")
    expect_true(validObject(out))
    expect_equal(nrow(chromData(out)), nrow(peak_tbl))

    ## Chromatographic data should have been restricted to the same identifiers
    out_cd <- chromData(out)
    expect_true(all(out_cd$dataOrigin %in% unique(peak_tbl$dataOrigin)))
    expect_true(all(out_cd$chromIndex %in% unique(peak_tbl$chromIndex)))

    ## rtMin/rtMax should match peak.table values
    expect_equal(out_cd$rtMin, peak_tbl$rtMin)
    expect_equal(out_cd$rtMax, peak_tbl$rtMax)

    bad_tbl <- data.frame(
        rtMin = 1, rtMax = 2,
        msLevel = 9L, dataOrigin = "nonexistent_file.mzML"
    )
    expect_error(
        chromExtract(be, bad_tbl, by = c("dataOrigin", "msLevel")),
        regexp = "do not exist|must be present"
    )
    pd <- peaksData(out) # should work
    expect_equal(length(pd), nrow(peak_tbl))
    expect_true(all(sapply(pd, nrow) > 0))

    expect_equal(length(.peaksData(out)), nrow(peak_tbl)) ## i think this fails

})

