test_that("ChromBackendSpectra works", {
    expect_true(isReadOnly(be_sp))
    expect_false(be_sp@inMemory)
    expect_false(identical(peaksData(be_sp), be_sp@peaksData))
    expect_true(identical(length(peaksData(be_sp)), length((be_sp@peaksData))))
    expect_true(identical(
        chromData(be_sp),
        fillCoreChromVariables(be_sp@chromData)
    ))
    expect_false(supportsSetBackend(be_sp))
})

test_that("backendInitialize works", {
    expect_equal(
        backendInitialize(ChromBackendSpectra(), spectra = Spectra()),
        ChromBackendSpectra()
    )
    expect_error(
        backendInitialize(ChromBackendSpectra(), spectra = numeric()),
        "must be a 'Spectra'"
    )
    expect_error(backendInitialize(ChromBackendSpectra(),
        spectra = s,
        factorize.by = c("nope", "nope2")
    ), "variables must exist in")
    expect_error(
        backendInitialize(ChromBackendSpectra(),
            summarize.method = "nope"
        ),
        "'arg' should be one of "
    )
    expect_error(
        backendInitialize(ChromBackendSpectra(),
            spectra = s,
            chromData = matrix()
        ),
        "must be a 'data.frame'"
    )
    expect_error(
        backendInitialize(ChromBackendSpectra(),
            spectra = s,
            chromData = data.frame(
                a = c(1, 2),
                b = c("yes", "no")
            )
        ),
        "All 'factorize.by' variables"
    )
    expect_error(
        backendInitialize(ChromBackendSpectra(),
            summarize.method = "mean"
        ),
        "should be one of"
    )

    df <- data.frame(
        chromIndex = 1:3, msLevel = 1:3, mz = 1:3,
        dataOrigin = dataOrigin(s)[1]
    )
    bd_tmp <- backendInitialize(ChromBackendSpectra(),
        spectra = s,
        chromData = df,
        factorize.by = c("msLevel", "dataOrigin")
    )
    df$chromSpectraIndex <- chromSpectraIndex(bd_tmp)
    expect_identical(bd_tmp@chromData, df)
})

test_that("show method for ChromBackendSpectra works correctly", {
    expect_output(show(be_sp), "1651 spectra")
    expect_output(show(be_sp), "3 chromatograms")
    expect_output(show(be_sp), "chromIndex")
    expect_output(show(be_sp), "msLevel")
    expect_output(show(be_sp), "mz")
    tmp <- be_sp
    tmp@inMemory <- TRUE
    expect_output(show(tmp), "Peaks data is cached in memory")
})

test_that("replacement method works", {
    tmp <- be_sp
    peaksData(tmp)[[1]] <- peaksData(be_sp)[[1]] + 1
    expect_false(identical(peaksData(tmp)[[1]], peaksData(be_sp)[[1]]))
    expect_true(identical(peaksData(tmp)[[1]], peaksData(be_sp)[[1]] + 1))
    expect_true(tmp@inMemory)
    expect_false(be_sp@inMemory)

    cd <- chromData(tmp)
    cd$mz <- 1
    expect_true(!identical(chromData(tmp), cd))
    chromData(tmp) <- cd
    expect_equal(chromData(tmp), cd)
})

test_that("backendParallelFactor works", {
    expect_equal(levels(backendParallelFactor(be_sp)), unique(chromSpectraIndex(be_sp)))
})

test_that("error message work", {
    expect_error(
        peaksData(be_sp, columns = "notacolumn"),
        "undefined columns selected"
    )
})

test_that("factorize() works", {
    expect_error(
        factorize(be_sp, factorize.by = "nope"),
        "variables must be in chromData"
    )
    expect_error(
        factorize(be_sp, factorize.by = "chromIndex"),
        "variables must be in the"
    )
    tmp <- be_sp
    tmp$msLevel <- c(1L, 2L, 3L)
    idx_before <- chromSpectraIndex(tmp)
    tmp <- factorize(tmp)
    idx_after <- chromSpectraIndex(tmp)
    expect_false(identical(idx_before, idx_after))
    expect_identical(
        levels(chromSpectraIndex(tmp)),
        levels(tmp@spectra$chromSpectraIndex)
    )
})

test_that("chromSpectraIndex works", {
    expect_error(
        chromSpectraIndex(1),
        "object must be a"
    )
    expect_equal(be_sp@chromData$chromSpectraIndex, chromSpectraIndex(be_sp))
})
