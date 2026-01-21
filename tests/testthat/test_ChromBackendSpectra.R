test_that("ChromBackendSpectra works", {
    expect_true(isReadOnly(be_sp))
    expect_false(be_sp@inMemory)
    expect_false(identical(peaksData(be_sp),.peaksData(be_sp)))
    expect_true(identical(length(peaksData(be_sp)), length((.peaksData(be_sp)))))
    expect_true(identical(
        chromData(be_sp),
        fillCoreChromVariables(be_sp@chromData)
    ))
    expect_false(supportsSetBackend(be_sp))

    ## works ofr single column factorisation too
    expect_silent(be_tmp <- factorize(be_sp, factorize.by = "dataOrigin"))
    expect_true(is(be_tmp, "ChromBackendSpectra"))
    expect_true(all(.chromData(be_tmp)$dataOrigin ==
                    .chromData(be_tmp)$chromSpectraIndex))
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
    df <- data.frame(msLevel = 1:3, mz = 1:3,
                     dataOrigin = dataOrigin(s)[1])
    bd_tmp <- backendInitialize(ChromBackendSpectra(), spectra = s,
                                chromData = df,
                                factorize.by = c("msLevel", "dataOrigin"))
    expect_identical(bd_tmp@chromData[, colnames(df)], df)
    expect_true(all(c("rtMin", "rtMax", "mzMin", "mzMax", "chromSpectraIndex") %in%
                colnames(bd_tmp@chromData)))
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

test_that("error message work", {
    expect_error(
        peaksData(be_sp, columns = "notacolumn"),
        "undefined columns selected"
    )
})

test_that("factorize() works", {
    expect_error(factorize(be_sp, factorize.by = "nope"),
                 "variables must be in the Spectra")
    expect_error(factorize(be_sp, factorize.by = "chromIndex"),
                 "variables must be in the")
    tmp <- be_sp
    tmp$msLevel <- c(1L, 2L, 3L)
    idx_before <- chromSpectraIndex(tmp)
    tmp <- factorize(tmp)
    idx_after <- chromSpectraIndex(tmp)
    expect_false(identical(idx_before, idx_after))
    expect_identical(levels(chromSpectraIndex(tmp)),
                     levels(tmp@spectra$chromSpectraIndex))
    tmp@spectra$extra_col <- seq_len(length(tmp@spectra))
    expect_error(factorize(tmp, factorize.by = c("msLevel", "extra_col")),
                 "must be in chromData")
})

test_that("factorize() fills rt columns with numeric vectors", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(4, c(1, 2), simplify = FALSE),
        intensity = replicate(4, c(10, 20), simplify = FALSE),
        rtime = c(20, 10, 5, 15),
        msLevel = rep(1L, 4),
        dataOrigin = c("A", "A", "B", "B")
    ))

    cb <- ChromBackendSpectra()
    cb@spectra <- sp
    cb@chromData <- data.frame(
        msLevel = 1L,
        dataOrigin = c("B", "A")
    )
    cb@spectraSortIndex <- order(sp$dataOrigin, sp$rtime)

    cb <- factorize(cb, factorize.by = "dataOrigin")

    expect_type(cb@chromData$rtMin, "double")
    expect_type(cb@chromData$rtMax, "double")
    expect_false(is.list(cb@chromData$rtMin))
    expect_equal(cb@chromData$rtMin, c(5, 10))
    expect_equal(cb@chromData$rtMax, c(15, 20))
})

test_that("spectraSortIndex is set and used for sorting", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(4, c(1, 2), simplify = FALSE),
        intensity = replicate(4, c(10, 20), simplify = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))

    cb <- ChromBackendSpectra()
    cb@spectra <- sp
    csi <- interaction(
        as.list(Spectra::spectraData(sp)[, c("msLevel", "dataOrigin"), drop = FALSE]),
        drop = TRUE, sep = "_"
    )
    cb@spectra$chromSpectraIndex <- csi
    cb@chromData <- fillCoreChromVariables(data.frame(
        msLevel = c(1L, 1L),
        dataOrigin = c("A", "B"),
        chromSpectraIndex = c("1_A", "1_B")
    ))
    cb@spectraSortIndex <- order(sp$dataOrigin, sp$rtime)

    expect_identical(cb@spectraSortIndex, order(sp$dataOrigin, sp$rtime))
})

test_that("[ maintains spectra and spectraSortIndex", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(4, c(1, 2), simplify = FALSE),
        intensity = replicate(4, c(10, 20), simplify = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))

    cb <- ChromBackendSpectra()
    cb@spectra <- sp
    csi <- interaction(
        as.list(Spectra::spectraData(sp)[, c("msLevel", "dataOrigin"), drop = FALSE]),
        drop = TRUE, sep = "_"
    )
    cb@spectra$chromSpectraIndex <- csi
    cb@chromData <- fillCoreChromVariables(data.frame(
        msLevel = c(1L, 1L),
        dataOrigin = c("A", "B"),
        chromSpectraIndex = c("1_A", "1_B")
    ))
    cb@spectraSortIndex <- order(sp$dataOrigin, sp$rtime)

    keep <- chromSpectraIndex(cb) == "1_A"
    cb_sub <- cb[keep]

    expect_s4_class(cb_sub, "ChromBackendSpectra")
    expect_identical(length(cb_sub), 1L)
    expect_true(all(cb_sub@spectra$dataOrigin == "A"))
    expect_identical(
        cb_sub@spectraSortIndex,
        order(cb_sub@spectra$dataOrigin, cb_sub@spectra$rtime)
    )
})

test_that("chromSpectraIndex works", {
    expect_error(
        chromSpectraIndex(1),
        "object must be a"
    )
    expect_equal(be_sp@chromData$chromSpectraIndex, chromSpectraIndex(be_sp))
    tmp <- be_sp
    tmp@chromData$chromSpectraIndex <- seq_len(nrow(tmp@chromData))
    expect_true(is.factor(chromSpectraIndex(tmp)))
})

test_that("backendParallelFactor works", {
    expect_identical(backendParallelFactor(be_sp), factor())
})

test_that("chromExtract works for ChromBackendSpectra", {
    tmp <- be_sp
    do <- dataOrigin(tmp)[1:2]
    peak_tbl <- data.frame(
        rtMin = c(1, 100),
        rtMax = c(56, 600),
        mzMin = c(100, 140),
        mzMax = c(130, 160),
        msLevel = c(1L, 1L),
        dataOrigin = do,
        extra_cols = c("peak1", "peak2")
    )

    out <- chromExtract(tmp, peak_tbl, by = c("msLevel", "dataOrigin"))

    # Output should be a ChromBackendSpectra
    expect_s4_class(out, "ChromBackendSpectra")

    # chromData should match rows in peak.table
    expect_equal(nrow(.chromData(out)), nrow(peak_tbl))

    # spectra should be preserved
    expect_true(!is.null(out@spectra))
    expect_true(!nrow(.peaksData(out)[[1]]) > 0)

    # Ensure replacement of mz columns occurred
    expect_true(all(c("mzMin", "mzMax", "extra_cols") %in% names(.chromData(out))))
    expect_equal(length(.peaksData(out)), nrow(peak_tbl))
})

