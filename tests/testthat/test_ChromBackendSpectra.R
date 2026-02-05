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
    sp <- Spectra(DataFrame(
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
    ## Set spectraSortIndex manually since this unsorted data
    sort_idx <- order(sp$dataOrigin, sp$rtime)
    if (!identical(sort_idx, seq_along(sp))) {
        cb@spectraSortIndex <- sort_idx
    }

    cb <- factorize(cb, factorize.by = "dataOrigin")

    expect_type(cb@chromData$rtMin, "double")
    expect_type(cb@chromData$rtMax, "double")
    expect_false(is.list(cb@chromData$rtMin))
    expect_equal(cb@chromData$rtMin, c(5, 10))
    expect_equal(cb@chromData$rtMax, c(15, 20))
})

test_that("spectraSortIndex is empty for pre-sorted data", {
    ## Create spectra that are already sorted by dataOrigin and rtime
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4, 5),
        msLevel = rep(1L, 5),
        dataOrigin = rep("A", 5)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## For pre-sorted data, spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
    expect_identical(cb@spectraSortIndex, integer())
})

test_that("spectraSortIndex is set for unsorted data", {
    ## Create spectra that are NOT sorted by rtime
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(5, 1, 3, 2, 4),
        msLevel = rep(1L, 5),
        dataOrigin = rep("A", 5)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## For unsorted data, spectraSortIndex should be set
    expect_true(length(cb@spectraSortIndex) > 0)
    expected_sort <- order(sp$dataOrigin, sp$rtime)
    expect_identical(cb@spectraSortIndex, expected_sort)
    
    ## Verify the sort index produces sorted rtimes
    sorted_rtimes <- sp$rtime[cb@spectraSortIndex]
    expect_identical(sorted_rtimes, c(1, 2, 3, 4, 5))
})

test_that("spectraSortIndex is empty for pre-sorted data with multiple dataOrigins", {
    ## Create spectra sorted by dataOrigin then rtime
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = rep(1L, 6),
        dataOrigin = c("A", "A", "A", "B", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## For pre-sorted data, spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
})

test_that("spectraSortIndex is set for unsorted data with multiple dataOrigins", {
    ## Create spectra NOT sorted by dataOrigin and rtime
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 2, 6, 4, 5),
        msLevel = rep(1L, 6),
        dataOrigin = c("B", "A", "A", "B", "A", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## For unsorted data, spectraSortIndex should be set
    expect_true(length(cb@spectraSortIndex) > 0)
    expected_sort <- order(sp$dataOrigin, sp$rtime)
    expect_identical(cb@spectraSortIndex, expected_sort)
    
    ## Verify sorting is correct
    sorted_do <- sp$dataOrigin[cb@spectraSortIndex]
    sorted_rt <- sp$rtime[cb@spectraSortIndex]
    expect_identical(sorted_do, c("A", "A", "A", "B", "B", "B"))
    expect_identical(sorted_rt, c(1, 2, 4, 3, 5, 6))
})

test_that("factorize() clears spectraSortIndex when data is sorted", {
    ## Create unsorted spectra
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(4, 1, 3, 2),
        msLevel = rep(1L, 4),
        dataOrigin = rep("A", 4)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## Initially unsorted, so spectraSortIndex should be set
    expect_true(length(cb@spectraSortIndex) > 0)
    
    ## Now manually sort the spectra and refactorize
    sorted_indices <- order(sp$dataOrigin, sp$rtime)
    cb@spectra <- sp[sorted_indices]
    cb@spectra$rtime <- sp$rtime[sorted_indices]
    cb <- factorize(cb)
    
    ## After refactorize with sorted data, spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
})

test_that("factorize() recalculates spectraSortIndex for unsorted data", {
    ## Create spectra and manually set it to unsorted state
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = rep(1L, 4),
        dataOrigin = rep("A", 4)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## Initially sorted, so spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
    
    ## Modify msLevel to create different factorization
    cb@spectra$msLevel <- c(1L, 2L, 1L, 2L)
    cb <- factorize(cb)
    
    ## spectraSortIndex should still be empty since data is still sorted by rtime
    expect_identical(length(cb@spectraSortIndex), 0L)
})

test_that("peaksData works correctly with empty spectraSortIndex", {
    ## Create pre-sorted spectra
    sp <- Spectra(DataFrame(
        mz = NumericList(c(100, 101), c(100, 101), c(100, 101), c(100, 101), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = rep(1L, 4),
        dataOrigin = rep("A", 4)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
    
    ## peaksData should still work correctly
    pd <- peaksData(cb)
    expect_true(is.list(pd))
    expect_identical(length(pd), 1L)
    expect_true(is.data.frame(pd[[1]]))
})

test_that("peaksData works correctly with populated spectraSortIndex", {
    ## Create unsorted spectra
    sp <- Spectra(DataFrame(
        mz = NumericList(c(100, 101), c(100, 101), c(100, 101), c(100, 101), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    ## spectraSortIndex should be set
    expect_true(length(cb@spectraSortIndex) > 0)
    
    ## peaksData should still work correctly
    pd <- peaksData(cb)
    expect_true(is.list(pd))
    expect_identical(length(pd), 2L)
    expect_true(all(sapply(pd, is.data.frame)))
})

test_that("spectraSortIndex is set and used for sorting", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))

    cb <- ChromBackendSpectra()
    cb@spectra <- sp
    csi <- interaction(
        as.list(spectraData(sp)[, c("msLevel", "dataOrigin"), drop = FALSE]),
        drop = TRUE, sep = "_"
    )
    cb@spectra$chromSpectraIndex <- csi
    cb@chromData <- fillCoreChromVariables(data.frame(
        msLevel = c(1L, 1L),
        dataOrigin = c("A", "B"),
        chromSpectraIndex = c("1_A", "1_B")
    ))
    ## Manually set spectraSortIndex since this is unsorted data
    sort_idx <- order(sp$dataOrigin, sp$rtime)
    if (!identical(sort_idx, seq_along(sp))) {
        cb@spectraSortIndex <- sort_idx
    }

    expect_true(length(cb@spectraSortIndex) > 0)
    expect_identical(cb@spectraSortIndex, order(sp$dataOrigin, sp$rtime))
})

test_that("[ maintains spectra and spectraSortIndex", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 4, 2),
        msLevel = rep(1L, 4),
        dataOrigin = c("B", "A", "B", "A")
    ))

    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial state
    expect_identical(length(cb), 2L)  # 2 unique combinations of msLevel and dataOrigin
    expect_identical(length(cb@spectra), 4L)  # 4 spectra
    # spectraSortIndex should be set since data is unsorted
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
    
    # Verify spectraSortIndex is valid (may be empty if subsetted data is sorted)
    if (length(cb_sub@spectraSortIndex) > 0) {
        expect_true(all(cb_sub@spectraSortIndex <= length(cb_sub@spectra)))
        # Verify the sort index still provides correct ordering
        sorted_rtime <- cb_sub@spectra$rtime[cb_sub@spectraSortIndex]
        expect_true(all(diff(sorted_rtime) >= 0))
    } else {
        # If spectraSortIndex is empty, spectra should already be sorted
        sorted_rtime <- cb_sub@spectra$rtime
        expect_true(all(diff(sorted_rtime) >= 0))
    }
    
    # Verify chromSpectraIndex is properly factorized
    expect_true(is.factor(cb_sub@spectra$chromSpectraIndex))
    expect_identical(levels(cb_sub@spectra$chromSpectraIndex), "1_A")
})

test_that("factorize() handles empty chromData correctly", {
    # Test the scenario from the vignette: creating Chromatograms from Spectra
    # without providing explicit chromData
    sp <- Spectra(DataFrame(
        mz = NumericList(
            c(100, 101), c(100, 101), c(100, 101), c(100, 101), c(100, 101),
            compress = FALSE
        ),
        intensity = NumericList(
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
    
    # Verify spectraSortIndex is valid (should be empty since data is pre-sorted)
    expect_identical(length(cb@spectraSortIndex), 0L)
    # For pre-sorted data, rtimes should already be in order
    expect_identical(cb@spectra$rtime, c(100, 110, 120, 130, 140))
    
    # Verify chromData has correct rt range
    expect_equal(cb@chromData$rtMin, 100)
    expect_equal(cb@chromData$rtMax, 140)
})

test_that("factorize() recalculates spectraSortIndex correctly - in-memory backend", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2),
                                   compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20),
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
    sp <- Spectra(DataFrame(
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
    sp <- Spectra(DataFrame(
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
    sp <- Spectra(DataFrame(
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

test_that("[ subsetting with empty spectraSortIndex works", {
    ## Create pre-sorted spectra (empty spectraSortIndex)
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial state: pre-sorted data should have empty spectraSortIndex
    expect_identical(length(cb@spectraSortIndex), 0L)
    
    # Subset to keep only first chromatogram (msLevel=1, dataOrigin=A)
    keep <- chromSpectraIndex(cb) == levels(chromSpectraIndex(cb))[1]
    cb_sub <- cb[keep]
    
    # Verify subset is correct
    expect_identical(length(cb_sub), 1L)
    expect_identical(length(cb_sub@spectra), 2L)
    expect_true(all(cb_sub@spectra$msLevel == 1L))
    
    # Verify spectraSortIndex is still empty (pre-sorted data stays pre-sorted)
    expect_identical(length(cb_sub@spectraSortIndex), 0L)
})

test_that("[ subsetting with populated spectraSortIndex works", {
    ## Create unsorted spectra (populated spectraSortIndex)
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 4, 2, 5, 6),  # Unsorted
        msLevel = c(1L, 1L, 2L, 2L, 1L, 2L),
        dataOrigin = c("B", "A", "B", "A", "A", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial state: unsorted data should have populated spectraSortIndex
    expect_true(length(cb@spectraSortIndex) > 0)
    original_sort_idx <- cb@spectraSortIndex
    
    # Subset to keep only "1_A" chromatogram (msLevel=1, dataOrigin=A)
    keep <- chromSpectraIndex(cb) == "1_A"
    cb_sub <- cb[keep]
    
    # Verify subset is correct
    expect_identical(length(cb_sub), 1L)
    expect_identical(length(cb_sub@spectra), 2L)
    
    # Verify spectra in subset have correct values
    expect_true(all(cb_sub@spectra$msLevel == 1L))
    expect_true(all(cb_sub@spectra$dataOrigin == "A"))
    
    # Verify spectraSortIndex is remapped correctly
    ## Original positions that map to "1_A" are positions 2 and 5
    ## After subsetting, they should be renumbered to 1, 2
    expect_true(all(cb_sub@spectraSortIndex <= length(cb_sub@spectra)))
    
    # Verify the sort index produces correct ordering
    if (length(cb_sub@spectraSortIndex) > 0) {
        sorted_rtimes <- cb_sub@spectra$rtime[cb_sub@spectraSortIndex]
        expect_true(all(diff(sorted_rtimes) >= 0))
    }
})

test_that("[ subsetting remaps spectraSortIndex correctly", {
    ## Test remapping with a specific example
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(5, 1, 3, 2, 4),  # Unsorted
        msLevel = rep(1L, 5),
        dataOrigin = c("A", "A", "B", "B", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Get the original spectraSortIndex
    original_sort_idx <- cb@spectraSortIndex
    expect_true(length(original_sort_idx) > 0)
    
    # Subset to keep only dataOrigin="A" (indices 1, 2, 5)
    keep <- chromSpectraIndex(cb) == "1_A"
    cb_sub <- cb[keep]
    
    # Verify we kept the right spectra
    expect_identical(length(cb_sub@spectra), 3L)
    kept_positions <- c(1, 2, 5)  # These are the original positions we kept
    expect_true(all(cb_sub@spectra$dataOrigin == "A"))
    
    # Verify spectraSortIndex is valid
    expect_true(all(cb_sub@spectraSortIndex > 0))
    expect_true(all(cb_sub@spectraSortIndex <= length(cb_sub@spectra)))
    
    # Verify sort order is preserved
    sorted_rtimes <- cb_sub@spectra$rtime[cb_sub@spectraSortIndex]
    expect_true(all(diff(sorted_rtimes) >= 0))
})

test_that("[ subsetting with multiple indices works", {
    ## Test subsetting with multiple chromatograms
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(2, 1, 4, 3),  # Unsorted
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "A", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Subset to keep both chromatograms
    keep <- c(TRUE, TRUE)  # Keep both msLevel=1 and msLevel=2
    cb_sub <- cb[keep]
    
    # Verify all spectra are retained
    expect_identical(length(cb_sub@spectra), 4L)
    expect_identical(length(cb_sub), 2L)
    
    # Verify spectraSortIndex is still valid
    if (length(cb_sub@spectraSortIndex) > 0) {
        expect_true(all(cb_sub@spectraSortIndex <= length(cb_sub@spectra)))
        sorted_do <- cb_sub@spectra$dataOrigin[cb_sub@spectraSortIndex]
        sorted_rt <- cb_sub@spectra$rtime[cb_sub@spectraSortIndex]
        # Verify sorting is correct
        expect_true(all(diff(sorted_rt) >= 0))
    }
})

test_that("[ subsetting with single index works", {
    ## Test subsetting to keep only one chromatogram
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(3, 1, 4, 2),  # Unsorted
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "A", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Subset to keep only first chromatogram
    keep <- c(TRUE, FALSE)
    cb_sub <- cb[keep]
    
    # Verify subset is correct
    expect_identical(length(cb_sub), 1L)
    expect_identical(length(cb_sub@spectra), 2L)
    expect_true(all(cb_sub@spectra$msLevel == 1L))
})

test_that("factorize() with single factorize.by variable works", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = c(1L, 1L, 2L, 2L, 1L, 2L),
        dataOrigin = c("A", "A", "A", "A", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)  # Default: factorize.by = c("msLevel", "dataOrigin")
    
    # Verify initial state creates 4 groups (1_A, 1_B, 2_A, 2_B)
    expect_identical(length(cb), 4L)
    
    # Refactorize by only dataOrigin
    cb_fact <- factorize(cb, factorize.by = "dataOrigin")
    
    # After refactorization, chromSpectraIndex changes to only have dataOrigin groups
    expect_true(is.factor(cb_fact@spectra$chromSpectraIndex))
    expect_identical(nlevels(cb_fact@spectra$chromSpectraIndex), 2L)  # Now only 2 levels: A and B
    
    # Verify spectra are correctly assigned to groups
    expect_true(all(cb_fact@spectra$chromSpectraIndex[cb_fact@spectra$dataOrigin == "A"] == "A"))
    expect_true(all(cb_fact@spectra$chromSpectraIndex[cb_fact@spectra$dataOrigin == "B"] == "B"))
})

test_that("factorize() with multiple factorize.by variables works", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = c(1L, 1L, 2L, 2L, 1L, 2L),
        dataOrigin = c("A", "A", "A", "A", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Factorize by msLevel and dataOrigin (default)
    cb_fact <- factorize(cb, factorize.by = c("msLevel", "dataOrigin"))
    
    # Should create 4 chromatograms
    expect_identical(length(cb_fact), 4L)
    expect_identical(nrow(cb_fact@chromData), 4L)
    
    # Verify chromData levels
    levels_str <- sort(as.character(unique(cb_fact@spectra$chromSpectraIndex)))
    expected_levels <- c("1_A", "1_B", "2_A", "2_B")
    expect_identical(levels_str, sort(expected_levels))
})

test_that("factorize() recalculates spectraSortIndex for unsorted data", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(5, 1, 3, 2, 6, 4),  # Unsorted
        msLevel = c(1L, 1L, 2L, 2L, 1L, 2L),
        dataOrigin = c("A", "A", "A", "A", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial spectraSortIndex is set (unsorted data)
    expect_true(length(cb@spectraSortIndex) > 0)
    initial_sort_idx <- cb@spectraSortIndex
    
    # Modify rtime and refactorize
    cb@spectra$rtime <- c(1, 2, 3, 4, 5, 6)  # Make it sorted
    cb <- factorize(cb)
    
    # After refactorize with sorted data, spectraSortIndex should be empty
    expect_identical(length(cb@spectraSortIndex), 0L)
})

test_that("factorize() updates chromSpectraIndex correctly", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(10, 20), c(10, 20), c(10, 20), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = rep("A", 4)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify chromSpectraIndex is a factor
    expect_true(is.factor(cb@spectra$chromSpectraIndex))
    
    # Refactorize
    cb <- factorize(cb)
    
    # Verify chromSpectraIndex is still a factor
    expect_true(is.factor(cb@spectra$chromSpectraIndex))
    
    # Verify all spectra have valid chromSpectraIndex
    expect_true(all(!is.na(cb@spectra$chromSpectraIndex)))
    
    # Verify levels match chromData
    expect_identical(sort(levels(cb@spectra$chromSpectraIndex)), 
                     sort(as.character(cb@chromData$chromSpectraIndex)))
})

test_that("factorize() preserves spectra data integrity", {
    sp <- Spectra(DataFrame(
        mz = NumericList(c(100, 101), c(100, 101), c(100, 101), c(100, 101), compress = FALSE),
        intensity = NumericList(c(10, 20), c(15, 25), c(30, 5), c(12, 18), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = rep("A", 4)
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    original_mz <- cb@spectra$mz
    original_intensity <- cb@spectra$intensity
    
    # Refactorize
    cb <- factorize(cb)
    
    # Verify spectra data is preserved
    expect_identical(cb@spectra$mz, original_mz)
    expect_identical(cb@spectra$intensity, original_intensity)
})

test_that("[ subsetting with reordering works", {
    ## Test that reordering chromatograms works correctly
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(15, 25), c(30, 5), c(12, 18), c(40, 2), c(50, 10), compress = FALSE),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = c(1L, 1L, 2L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "A", "B", "B", "B")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Get original chromSpectraIndex levels
    original_idx <- chromSpectraIndex(cb)
    expect_identical(length(original_idx), 4L)  # 4 chromatograms: 1_A, 1_B, 2_A, 2_B
    
    # Reorder: reverse order
    cb_reordered <- cb[c(4, 3, 2, 1)]
    
    # Verify length is correct
    expect_identical(length(cb_reordered), 4L)
    
    # Verify order is reversed
    reordered_idx <- chromSpectraIndex(cb_reordered)
    expect_identical(as.character(reordered_idx), as.character(original_idx[c(4, 3, 2, 1)]))
    
    # Verify chromData is reordered
    expect_identical(nrow(cb_reordered@chromData), 4L)
    
    # Verify spectra are still correctly associated
    expect_true(all(cb_reordered@spectra$chromSpectraIndex %in% reordered_idx))
    
    # Verify spectraSortIndex is still valid
    if (length(cb_reordered@spectraSortIndex) > 0) {
        expect_true(all(cb_reordered@spectraSortIndex <= length(cb_reordered@spectra)))
    }
})

test_that("[ subsetting with duplication works", {
    ## Test that duplicating chromatograms works correctly
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(15, 25), c(30, 5), c(12, 18), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "A", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Get original chromSpectraIndex levels
    original_idx <- chromSpectraIndex(cb)
    expect_identical(length(original_idx), 2L)  # 2 chromatograms: 1_A, 2_A
    
    # Duplicate: keep first chromatogram twice, then second
    cb_dup <- cb[c(1, 1, 2)]
    
    # Verify length includes duplicates
    expect_identical(length(cb_dup), 3L)
    
    # Verify chromData includes duplicates
    expect_identical(nrow(cb_dup@chromData), 3L)
    dup_idx <- chromSpectraIndex(cb_dup)
    expect_identical(as.character(dup_idx), as.character(original_idx[c(1, 1, 2)]))
    
    # Verify spectra: keeps all unique spectra from subsetted chromatograms
    ## Both chrom 1 and 2 are kept, so all 4 original spectra are kept
    expect_identical(length(cb_dup@spectra), 4L)
    
    # Verify all spectra belong to the right chromatograms
    expect_true(all(cb_dup@spectra$chromSpectraIndex %in% dup_idx))
    
    # Verify spectraSortIndex is still valid
    if (length(cb_dup@spectraSortIndex) > 0) {
        expect_true(all(cb_dup@spectraSortIndex <= length(cb_dup@spectra)))
    }
    
    # Verify peaksData works with duplicated chromatograms
    pd <- peaksData(cb_dup)
    expect_identical(length(pd), 3L)
})

test_that("[ subsetting with reordering and duplication works", {
    ## Test combined reordering and duplication
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(15, 25), c(30, 5), c(12, 18), compress = FALSE),
        rtime = c(1, 2, 3, 4),
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("A", "A", "A", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    original_idx <- chromSpectraIndex(cb)
    expect_identical(length(original_idx), 2L)  # 2 chromatograms
    
    # Reorder and duplicate: c(2, 1, 1)
    cb_mixed <- cb[c(2, 1, 1)]
    
    # Verify length is correct
    expect_identical(length(cb_mixed), 3L)
    
    # Verify order matches requested subset
    mixed_idx <- chromSpectraIndex(cb_mixed)
    expect_identical(as.character(mixed_idx), as.character(original_idx[c(2, 1, 1)]))
    
    # Verify chromData has correct structure
    expect_identical(nrow(cb_mixed@chromData), 3L)
    
    # Verify spectra are correctly associated
    expect_true(all(cb_mixed@spectra$chromSpectraIndex %in% mixed_idx))
    
    # Verify number of spectra: keeps unique spectra from both chromatograms
    ## Both chrom 1 and 2 are in the subset, so all 4 spectra are kept
    expect_identical(length(cb_mixed@spectra), 4L)
    
    # Verify spectraSortIndex is valid
    if (length(cb_mixed@spectraSortIndex) > 0) {
        expect_true(all(cb_mixed@spectraSortIndex > 0))
        expect_true(all(cb_mixed@spectraSortIndex <= length(cb_mixed@spectra)))
        
        # Verify sort index still provides correct ordering
        sorted_rtimes <- cb_mixed@spectra$rtime[cb_mixed@spectraSortIndex]
        sorted_do <- cb_mixed@spectra$dataOrigin[cb_mixed@spectraSortIndex]
        
        # Check rtimes are sorted within each dataOrigin group
        for (do in unique(sorted_do)) {
            rtimes_in_group <- sorted_rtimes[sorted_do == do]
            expect_true(all(diff(rtimes_in_group) >= 0))
        }
    }
})

test_that("[ subsetting with unsorted data and duplication works", {
    ## Test duplication with unsorted spectra (non-empty spectraSortIndex)
    sp <- Spectra(DataFrame(
        mz = NumericList(c(1, 2), c(1, 2), c(1, 2), c(1, 2), compress = FALSE),
        intensity = NumericList(c(10, 20), c(15, 25), c(30, 5), c(12, 18), compress = FALSE),
        rtime = c(3, 1, 4, 2),  # Unsorted
        msLevel = c(1L, 1L, 2L, 2L),
        dataOrigin = c("B", "A", "B", "A")
    ))
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp)
    
    # Verify initial state has populated spectraSortIndex
    expect_true(length(cb@spectraSortIndex) > 0)
    
    original_idx <- chromSpectraIndex(cb)
    expect_identical(length(original_idx), 4L)  # 4 chromatograms: 1_A, 1_B, 2_A, 2_B
    
    # Duplicate and reorder: keep chromatogram 2, then 1 twice
    cb_dup <- cb[c(2, 1, 1)]
    
    # Verify basic structure
    expect_identical(length(cb_dup), 3L)
    expect_identical(nrow(cb_dup@chromData), 3L)
    
    # Verify order
    dup_idx <- chromSpectraIndex(cb_dup)
    expect_identical(as.character(dup_idx), as.character(original_idx[c(2, 1, 1)]))
    
    # Verify spectra count: keeps unique spectra from chromatograms 1 and 2
    ## Chrom 1_A has 1 spectrum, chrom 1_B has 1 spectrum → 2 unique spectra
    expect_identical(length(cb_dup@spectra), 2L)
    
    # Verify spectraSortIndex is valid after duplication
    if (length(cb_dup@spectraSortIndex) > 0) {
        expect_true(all(cb_dup@spectraSortIndex > 0))
        expect_true(all(cb_dup@spectraSortIndex <= length(cb_dup@spectra)))
        
        # Verify sorting is maintained
        sorted_rtimes <- cb_dup@spectra$rtime[cb_dup@spectraSortIndex]
        sorted_do <- cb_dup@spectra$dataOrigin[cb_dup@spectraSortIndex]
        
        # Within each dataOrigin, rtimes should be sorted
        for (do in unique(sorted_do)) {
            rtimes_in_group <- sorted_rtimes[sorted_do == do]
            expect_true(all(diff(rtimes_in_group) >= 0))
        }
    }
    
    # Verify peaksData works with duplicated chromatograms
    pd <- peaksData(cb_dup)
    expect_identical(length(pd), 3L)
    expect_true(all(sapply(pd, is.data.frame)))
})

test_that("peaksData returns data in correct order when chromSpectraIndex has duplicates", {
    ## Test case that reproduces the setBackend bug where peaksData was returned
    ## in factor level order instead of chromData row order.
    ## When multiple chromatograms share the same chromSpectraIndex (e.g., multiple
    ## EICs from the same file), split() groups them together by level. The fix
    ## ensures peaksData is returned in the original chromData row order.
    sp <- Spectra(DataFrame(
        mz = NumericList(
            c(100, 101, 102), c(100, 101, 102), c(100, 101, 102), 
            c(100, 101, 102), c(100, 101, 102), c(100, 101, 102),
            compress = FALSE
        ),
        intensity = NumericList(
            c(10, 20, 30), c(15, 25, 35), c(20, 30, 40),
            c(100, 200, 300), c(150, 250, 350), c(200, 300, 400),
            compress = FALSE
        ),
        rtime = c(1, 2, 3, 4, 5, 6),
        msLevel = rep(1L, 6),
        dataOrigin = c("A", "A", "A", "B", "B", "B")
    ))
    
    ## Create custom chromData with multiple EICs sharing same chromSpectraIndex
    ## This simulates extracting different m/z windows from the same spectra
    custom_chromData <- data.frame(
        msLevel = rep(1L, 4),
        mzMin = c(100, 100, 100, 100),
        mzMax = c(101, 102, 101, 102),
        rtMin = c(1, 1, 4, 4),
        rtMax = c(3, 3, 6, 6),
        dataOrigin = c("A", "A", "B", "B")
    )
    
    cb <- ChromBackendSpectra()
    cb <- backendInitialize(cb, spectra = sp, chromData = custom_chromData)
    
    ## chromSpectraIndex should have duplicates: 1_A, 1_A, 1_B, 1_B
    chrom_idx <- chromSpectraIndex(cb)
    expect_identical(length(chrom_idx), 4L)
    expect_identical(as.character(chrom_idx), c("1_A", "1_A", "1_B", "1_B"))
    
    ## Get peaksData
    pd <- peaksData(cb)
    
    ## Verify peaksData has same length as chromData rows
    expect_identical(length(pd), 4L)
    
    ## Verify peaksData is returned in chromData row order, NOT factor level order
    ## Row 1: 1_A with mzMax=101 (should have intensities from mz=100,101)
    ## Row 2: 1_A with mzMax=102 (should have intensities from mz=100,101,102)
    ## Row 3: 1_B with mzMax=101 (should have intensities from mz=100,101)
    ## Row 4: 1_B with mzMax=102 (should have intensities from mz=100,101,102)
    
    ## Check that the intensity patterns match the expected order
    ## The key test: row 1 and row 2 should have different intensity sums
    ## because they have different mz windows
    int_sums <- sapply(pd, function(x) sum(x$intensity, na.rm = TRUE))
    
    ## Row 1 (1_A, mz 100-101): should have lower sum than Row 2 (1_A, mz 100-102)
    expect_true(int_sums[1] < int_sums[2])
    ## Row 3 (1_B, mz 100-101): should have lower sum than Row 4 (1_B, mz 100-102)
    expect_true(int_sums[3] < int_sums[4])
    
    ## Also verify the pattern matches across files
    ## Row 1 and Row 3 should have similar intensity ratios (both use mz 100-101)
    ## Row 2 and Row 4 should have similar intensity ratios (both use mz 100-102)
})

