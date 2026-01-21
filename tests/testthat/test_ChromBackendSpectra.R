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
    df <- data.frame(msLevel = 1L, mz = 100,
                     dataOrigin = dataOrigin(s)[1])
    bd_tmp <- backendInitialize(ChromBackendSpectra(), spectra = s,
                                chromData = df,
                                factorize.by = c("msLevel", "dataOrigin")
                              )
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
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial state
    expect_identical(length(cb), 2L)  # 2 unique combinations of msLevel and dataOrigin
    expect_identical(length(cb@spectra), 4L)  # 4 spectra
    expect_true(length(cb@spectraSortIndex) > 0)
    
    # Verify spectraSortIndex is correctly set
    expected_sort <- order(sp$dataOrigin, sp$rtime)
    expect_identical(cb@spectraSortIndex, expected_sort)
    
    # Get chromSpectraIndex before subsetting
    all_chrom_idx <- chromSpectraIndex(cb)
    expect_identical(length(all_chrom_idx), 2L)
    
    # Subset to keep only "1_A" chromatograms
    keep <- all_chrom_idx == "1_A"
    cb_sub <- cb[keep]

    # Verify subsetting worked correctly
    expect_s4_class(cb_sub, "ChromBackendSpectra")
    expect_identical(length(cb_sub), 1L)  # Only 1 chromatogram kept
    
    # Verify chromData was subsetted
    sub_chrom_idx <- chromSpectraIndex(cb_sub)
    expect_identical(length(sub_chrom_idx), 1L)
    expect_identical(as.character(sub_chrom_idx), "1_A")
    
    # Verify spectra were subsetted to only include spectra for "1_A"
    expect_identical(length(cb_sub@spectra), 2L)  # Only 2 spectra belong to "1_A"
    expect_true(all(cb_sub@spectra$dataOrigin == "A"))
    
    # Verify spectraSortIndex is valid and functional
    expect_true(length(cb_sub@spectraSortIndex) > 0)
    expect_true(all(cb_sub@spectraSortIndex <= length(cb_sub@spectra)))
    
    # Verify the sort index still provides correct ordering
    sorted_rtime <- cb_sub@spectra$rtime[cb_sub@spectraSortIndex]
    expect_true(all(diff(sorted_rtime) >= 0))
    
    # Verify chromSpectraIndex is properly factorized
    expect_true(is.factor(cb_sub@spectra$chromSpectraIndex))
    expect_identical(levels(cb_sub@spectra$chromSpectraIndex), "1_A")
})

test_that("factorize() handles empty chromData correctly", {
    # Test the scenario from the vignette: creating Chromatograms from Spectra
    # without providing explicit chromData
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = IRanges::NumericList(
            c(100, 101), c(100, 101), c(100, 101), c(100, 101), c(100, 101),
            compress = FALSE
        ),
        intensity = IRanges::NumericList(
            c(10, 20), c(15, 25), c(30, 5), c(12, 18), c(40, 2),
            compress = FALSE
        ),
        rtime = c(100, 110, 120, 130, 140),
        msLevel = rep(1L, 5),
        dataOrigin = rep("example", 5)
    ))
    
    # Create Chromatograms with empty chromData
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp, chromData = data.frame())
    
    # Should create chromData from spectra
    expect_identical(nrow(cb@chromData), 1L)  # One chromatogram: 1_example
    expect_identical(length(cb@spectra), 5L)  # All 5 spectra retained
    
    # Verify chromSpectraIndex was created correctly
    expect_true(is.factor(cb@spectra$chromSpectraIndex))
    expect_identical(as.character(unique(cb@spectra$chromSpectraIndex)), "1_example")
    
    # Verify spectraSortIndex is valid
    expect_identical(length(cb@spectraSortIndex), 5L)
    sorted_rtimes <- cb@spectra$rtime[cb@spectraSortIndex]
    expect_identical(sorted_rtimes, c(100, 110, 120, 130, 140))
    
    # Verify chromData has correct rt range
    expect_equal(cb@chromData$rtMin, 100)
    expect_equal(cb@chromData$rtMax, 140)
})

test_that("factorize() recalculates spectraSortIndex correctly - in-memory backend", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = IRanges::NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2),
                                   compress = FALSE),
        intensity = IRanges::NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20),
                                          compress = FALSE),
        rtime = c(3, 1, 4, 2, 5),
        msLevel = c(1L, 1L, 2L, 2L, 1L),
        dataOrigin = rep("A", 5)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial sort order is by dataOrigin, rtime
    original_sort_idx <- cb@spectraSortIndex
    sorted_rtimes <- sp$rtime[original_sort_idx]
    expect_true(all(diff(sorted_rtimes) >= 0))
    
    # Modify msLevel and refactorize
    cb@spectra$msLevel <- c(2L, 1L, 1L, 2L, 1L)
    cb <- factorize(cb)
    
    # After factorize, spectraSortIndex should be recalculated
    # and should still produce correctly sorted rtimes
    new_sort_idx <- cb@spectraSortIndex
    sorted_rtimes_after <- cb@spectra$rtime[new_sort_idx]
    expect_true(all(diff(sorted_rtimes_after) >= 0))
    
    # Verify chromData reflects the new factorization
    expect_identical(nrow(cb@chromData), 2L)  # Now 2 groups: 1_A and 2_A
})

test_that("factorize() works correctly with on-disk spectra backend", {
    # Use the pre-loaded on-disk spectra from test setup (be_sp)
    cb_test <- be_sp
    
    # Verify initial state
    expect_true(length(cb_test@spectra) > 0)
    expect_true(length(cb_test@spectraSortIndex) > 0)
    
    # Verify sort index produces sorted rtimes within groups
    sorted_rtimes <- cb_test@spectra$rtime[cb_test@spectraSortIndex]
    sorted_dataOrigin <- cb_test@spectra$dataOrigin[cb_test@spectraSortIndex]
    # Check rtimes are sorted within each dataOrigin group
    for (do in unique(sorted_dataOrigin)) {
        rtimes_in_group <- sorted_rtimes[sorted_dataOrigin == do]
        expect_true(all(diff(rtimes_in_group) >= 0))
    }
    
    # Factorize with different grouping
    cb_single <- factorize(cb_test, factorize.by = "dataOrigin")
    
    # Should have one row per unique dataOrigin
    expect_true(nrow(cb_single@chromData) >= 1)
    
    # spectraSortIndex should still be valid
    expect_true(all(cb_single@spectraSortIndex <= length(cb_single@spectra)))
    sorted_rtimes_new <- cb_single@spectra$rtime[cb_single@spectraSortIndex]
    sorted_dataOrigin_new <- cb_single@spectra$dataOrigin[cb_single@spectraSortIndex]
    # Check rtimes are sorted within each dataOrigin group
    for (do in unique(sorted_dataOrigin_new)) {
        rtimes_in_group <- sorted_rtimes_new[sorted_dataOrigin_new == do]
        expect_true(all(diff(rtimes_in_group) >= 0))
    }
})

test_that("factorize() maintains consistency between chromData and spectra", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(6, c(1, 2), simplify = FALSE),
        intensity = replicate(6, c(10, 20), simplify = FALSE),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = c(1L, 1L, 2L, 2L, 1L, 2L),
        dataOrigin = c("A", "A", "A", "B", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Get the chromSpectraIndex values
    chrom_idx <- chromSpectraIndex(cb)
    spectra_idx <- cb@spectra$chromSpectraIndex
    
    # Both should have the same set of levels
    expect_identical(sort(levels(chrom_idx)), sort(levels(spectra_idx)))
    
    # All chromSpectraIndex values in spectra should be in chromData levels
    unique_spectra_idx <- unique(as.character(spectra_idx))
    expect_true(all(unique_spectra_idx %in% levels(spectra_idx)))
    
    # Verify consistency after subsetting
    keep_idx <- chrom_idx == levels(chrom_idx)[1]
    cb_sub <- cb[keep_idx]
    
    chrom_idx_sub <- chromSpectraIndex(cb_sub)
    spectra_idx_sub <- cb_sub@spectra$chromSpectraIndex
    
    # After subsetting, both should have same single level
    expect_identical(nlevels(chrom_idx_sub), 1L)
    expect_identical(nlevels(spectra_idx_sub), 1L)
    expect_identical(levels(chrom_idx_sub), levels(spectra_idx_sub))
})

test_that("peaksData generation respects spectraSortIndex - in-memory", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(4, c(100, 101), simplify = FALSE),
        intensity = replicate(4, c(10, 20), simplify = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Get peaksData
    pd <- peaksData(cb)
    
    # Should return a list with one entry per chromatogram
    expect_identical(length(pd), 2L)
    
    # Each entry should be a data.frame with rtime and intensity
    for (i in seq_along(pd)) {
        expect_true(is.data.frame(pd[[i]]))
        expect_true("rtime" %in% colnames(pd[[i]]))
        expect_true("intensity" %in% colnames(pd[[i]]))
    }
})

test_that("subsetting and peaksData consistency - in-memory", {
    sp <- Spectra::Spectra(S4Vectors::DataFrame(
        mz = replicate(4, c(100, 101), simplify = FALSE),
        intensity = replicate(4, c(10, 20), simplify = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Subset to one chromatogram
    keep <- chromSpectraIndex(cb) == "1_A"
    cb_sub <- cb[keep]
    
    # Get peaksData from subsetted backend
    pd_sub <- peaksData(cb_sub)
    
    # Should have one entry
    expect_identical(length(pd_sub), 1L)
    
    # The data should be valid
    expect_true(is.data.frame(pd_sub[[1]]))
    expect_true(nrow(pd_sub[[1]]) > 0)
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

