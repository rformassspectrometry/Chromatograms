test_that("ChromBackendMemory function works", {
    expect_true(is.function(ChromBackendMemory))
    expect_true(is(ChromBackendMemory(), "ChromBackendMemory"))
    expect_true(supportsSetBackend(ChromBackendMemory()))
})

test_that("backendInitialize, ChromBackendMemory works", {
    expect_false(isReadOnly(be_empty))
    expect_error(
        backendInitialize(be_empty, chromData = matrix()),
        "needs to be a"
    )
    fill_test <- backendInitialize(be_empty, chromData = data.frame())
    expect_equal(.chromData(fill_test), fillCoreChromVariables(data.frame()))

    ### empty backend
    expect_true(is(be_empty, "ChromBackendMemory"))
    expect_true(is(.chromData(be_empty), "data.frame"))
    expect_true(nrow(.chromData(be_empty)) == 0)
    expect_true(all(names(.chromData(be_empty)) ==
                        names(.CORE_CHROM_VARIABLES)))
    expect_true(is.null(validChromData(.chromData(be_empty))))
    expect_true(is(.peaksData(be_empty), "list"))
    expect_true(length(.peaksData(be_empty)) == 1)
    expect_true(all(names(.peaksData(be_empty)[[1]]) ==
                        names(.CORE_PEAKS_VARIABLES)))
    expect_true(is.null(unlist(validChromData(.peaksData(be_empty)))))

    ### empty peaksData
    expect_true(is(be_cd, "ChromBackendMemory"))
    expect_true(is(.chromData(be_cd), "data.frame"))
    expect_true(nrow(.chromData(be_cd)) == nrow(cdata))
    expect_true(all(names(.chromData(be_cd)) %in% c(names(cdata),
                                                    "dataOrigin")))
    expect_true(is.null(validChromData(.chromData(be_cd))))
    expect_true(all(be_cd$dataOrigin == c("mem1", "mem2", "mem3")))

    expect_true(is(.peaksData(be_cd), "list"))
    expect_true(length(.peaksData(be_cd)) == nrow(.chromData(be_cd)))
    expect_true(all(vapply(.peaksData(be_cd), nrow, integer(1)) == 0))
    expect_true(all(vapply(.peaksData(be_cd), is.data.frame, logical(1))))
    expect_true(all(names(.peaksData(be_cd)[[1]]) %in%
                        names(.CORE_PEAKS_VARIABLES)))
    expect_true(is.null(unlist(validPeaksData(.peaksData(be_cd)))))

    ### full backend
    expect_true(is(be, "ChromBackendMemory"))
    expect_true(is(.chromData(be), "data.frame"))
    expect_true(nrow(.chromData(be)) == nrow(cdata))
    expect_true(all(names(.chromData(be)) %in% c(names(cdata), "dataOrigin")))
    expect_true(is.null(validChromData(.chromData(be))))

    expect_true(is(.peaksData(be), "list"))
    expect_true(length(.peaksData(be)) == nrow(.chromData(be)))
    expect_true(all(vapply(.peaksData(be), nrow, integer(1)) == nrow(pdata)))
    expect_true(all(vapply(.peaksData(be), is.data.frame, logical(1))))
    expect_true(all(names(.peaksData(be)[[1]]) %in%
                        names(.CORE_PEAKS_VARIABLES)))
    expect_true(is.null(unlist(validPeaksData(.peaksData(be)))))
})

test_that("backend removes NA columns", {
    chromData <- chromData(be)
    chromData$NAcol <- NA
    be2 <- backendInitialize(be, chromData)
    expect_false("NAcol" %in% names(.chromData(be2)))
})

test_that("backendMerge, ChromBackendMemory works", {
    be_merge <- backendMerge(c(be_cd, be))
    expect_true(is(be_merge, "ChromBackendMemory"))
    expect_true(is(.chromData(be_merge), "data.frame"))
    expect_true(is(.peaksData(be_merge), "list"))
    expect_true(nrow(.chromData(be_merge)) ==
                    nrow(.chromData(be_cd)) + nrow(.chromData(be)))
    expect_true(length(.peaksData(be_merge))
                == length(.peaksData(be_cd)) + length(.peaksData(be)))
    expect_equal(backendMerge(c(be_empty, be_empty)), be_empty)
})

test_that("show,ChromBackendMemory works", {
    expect_output(show(be_empty), "ChromBackendMemory")
    expect_output(show(be_cd), "3 chromatograms")
    expect_output(show(be_cd), "chromIndex")
    expect_output(show(be_cd), "msLevel")
    expect_output(show(be_cd), "mz")
})

test_that("chromData", {
    expect_error(chromData(be, "not a column"), "variables are not")
    expect_error(chromData(be) <- data.frame(), " rows")
    expect_error(chromData(be) <- matrix(), "is expected to be a 'data.frame'")
})

test_that("peaksData", {
    expect_error(peaksData(be, "not a column"), "variables are not")
    expect_error(peaksData(be) <- list(), " elements")
    expect_error(peaksData(be) <- matrix(), "is expected to be a list")
})

test_that("$<-,ChromBackendMemory works", {
    expect_error(
        be_cd$intensity <- c(100, 200),
        "length of 'value' needs to match the number of chromatograms in object."
    )

    new_peaks_intensity <- list(
        c(100, 200, 300, 400),
        c(50, 60),
        c(150, 250, 350, 450)
    )
    be$intensity <- new_peaks_intensity
    expect_equal(.peaksData(be)[[1]]$intensity, new_peaks_intensity[[1]])
    expect_equal(.peaksData(be)[[2]]$intensity, new_peaks_intensity[[2]])
    expect_equal(.peaksData(be)[[3]]$intensity, new_peaks_intensity[[3]])

    expect_error(
        be$intensity <- c(100, 200, 300),
        "The value for peaksData should be a list"
    )

    be_cd$new_var <- c("A", "B", "C")
    expect_equal(.chromData(be_cd)$new_var, c("A", "B", "C"))

    be_cd$mz <- c(111.1, 222.2, 333.3)
    expect_equal(.chromData(be_cd)$mz, c(111.1, 222.2, 333.3))
})

test_that("filterChromData works", {
    tmp <- be_cd
    tmp$new_var <- c(1, 2, 3)
    res <- filterChromData(tmp,
        variables = c("mz", "new_var"),
        ranges = c(134, 150, 1, 1),
        match = "any",
        keep = FALSE
    )

    expect_equal(nrow(chromData(res)), 1)
    expect_equal(chromData(res)$mz, 123.3)
    expect_error(
        filterChromData(tmp,
            variables = c("mz"),
            ranges = c(134), # Wrong length
            match = "any",
            keep = TRUE
        ),
        "needs to be twice the length"
    )
    expect_error(
        filterChromData(tmp,
            variables = c("invalid_var"),
            ranges = c(134, 150, 1, 1),
            match = "any",
            keep = TRUE
        ),
        "One or more values passed"
    )
    res <- filterChromData(tmp,
        variables = c("mz"),
        ranges = c(134, 150),
        match = "all",
        keep = FALSE
    )

    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 123.3))

    expect_identical(
        filterChromData(tmp,
            match = "all",
            keep = TRUE
        ),
        tmp
    )
    expect_error(
        filterChromData(tmp,
            variables = c("mz"),
            ranges = c("a", "b"),
            keep = TRUE
        ),
        "only support filtering for numerical"
    )
    expect_error(
        filterChromData(tmp,
            variables = c(200),
            ranges = c(134, 150),
            keep = TRUE
        ),
        "parameter needs to be of type"
    )
    res <- filterChromData(tmp,
        variables = c("mz", "new_var"),
        ranges = c(134, 150, 1, 1),
        match = "any",
        keep = TRUE
    )
    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 134.4))
    expect_equal(chromData(res)$new_var, c(1, 3))

    res <- filterChromData(tmp,
        variables = c("mz"),
        ranges = c(999, 1000), # No matching range
        match = "any",
        keep = FALSE
    )

    expect_identical(res, tmp)
})


test_that("split,ChrombackendMemory works", {
    be_split <- be
    be_split$new_vars <- c("a", "b", "b")
    tmp <- be
    tmp$newvar2s <- c(1, 2, 3)
    f <- factor(be_split$new_vars)
    res <- split(tmp, f)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "ChromBackendMemory")
    expect_s4_class(res[[2L]], "ChromBackendMemory")
    expect_equal(res[[1L]]$newvar2s, c(1))
    expect_equal(res[[2L]]$newvar2s, c(2, 3))

    res <- split(tmp, factor(be_split$new_vars, levels = c("b", "a")))
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "ChromBackendMemory")
    expect_s4_class(res[[2L]], "ChromBackendMemory")
    expect_equal(res[[2L]]$newvar2s, c(1))
    expect_equal(res[[1L]]$newvar2s, c(2, 3))
})


test_that("chromExtract works correctly for ChromBackendMemory", {
    peak_tbl <- data.frame(
        rtMin = c(12.3, 45.0, 12.5),
        rtMax = c(13.5, 46.3, 14.8),
        msLevel = c(1L, 1L, 1L),
        dataOrigin = c("mem1", "mem2", "mem3"),
        extra_cols = c(1, 2, 3)
    )

    out <- chromExtract(be, peak_tbl, by = c("msLevel", "dataOrigin"))

    expect_s4_class(out, "ChromBackendMemory")
    expect_true(validObject(out))

    ## chromData should have same number of rows as peak.table
    expect_equal(nrow(chromData(out)), nrow(peak_tbl))

    ## peaksData should be filtered per retention window
    pd_out <- peaksData(out)
    expect_type(pd_out, "list")
    expect_length(pd_out, 3)

    ## Check retention time filtering on first chromatogram
    rt1 <- pdata[[1]]$rtime
    in_range <- rt1 >= peak_tbl$rtMin[1] & rt1 <= peak_tbl$rtMax[1]
    expect_equal(pd_out[[1]]$rtime, rt1[in_range])

    ## intensities align with filtered rt
    expect_equal(length(pd_out[[1]]$rtime), length(pd_out[[1]]$intensity))

    ## Ensure no chromatograms missing
    expect_false(any(vapply(pd_out, nrow, numeric(1)) == 0))

    ## edge case: no match
    bad_tbl <- data.frame(
        rtMin = 1, rtMax = 2,
        msLevel = 2L,
        dataOrigin = "nonexistent"
    )
    expect_warning_or_error <- function(expr) {
        expect_error(expr, regexp = "do not exist|must be present",
                     fixed = FALSE)
    }
    expect_warning_or_error(chromExtract(be, bad_tbl,
                                         by = c("msLevel", "dataOrigin")))

})

test_that("chromExtract handles empty peaksData gracefully", {
    cdata <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100.1, 110.2),
        dataOrigin = c("x", "y")
    )
    be <- backendInitialize(new("ChromBackendMemory"), chromData = cdata)

    peak_tbl <- data.frame(
        rtMin = c(1, 2),
        rtMax = c(3, 4),
        msLevel = c(1L, 1L),
        dataOrigin = c("x", "y")
    )

    out <- chromExtract(be, peak_tbl, by = c("msLevel", "dataOrigin"))
    expect_s4_class(out, "ChromBackendMemory")
    expect_equal(nrow(chromData(out)), 2L)
    expect_equal(length(peaksData(out)), 2L)
    expect_true(all(vapply(peaksData(out), nrow, numeric(1)) == 0))
})

test_that("imputePeaksData() works correctly for ChromBackendMemory", {
    pdata <- list(
        data.frame(
            rtime = seq(12, 20, by = 0.5),
            intensity = c(123.3, 153.6, NA, NA, NA, 230, NA, 250, 260, 280,
                          300, 320, 340, 360, 380, 400, 420)
        ),
        data.frame(
            rtime = seq(45, 55, by = 0.5),  # length 21
            intensity = c(100, NA, 120, 130, 140, NA, NA, 170, 180, 190, 200,
                          210, 220, 230, 240, 250, 260, 270, 280, 290, 300)
        ),
        data.frame(
            rtime = seq(10, 18, by = 0.5),  # length 17
            intensity = c(NA, 153.6, 2354.3, 243.4, NA, NA, 280, 300, 320,
                          340, 360, 380, 400, 420, 440, NA, NA)
        )
    )

    tmp <- backendInitialize(be_empty, chromData = cdata, peaksData = pdata)

    expect_true(any(is.na(unlist(intensity(tmp)))))

    ## check each method
    methods <- c("linear", "spline")
    for (m in methods) {
        expect_silent({
            be_imp <- imputePeaksData(tmp, method = m)
        })

        ## No more NAs
        expect_false(any(is.na(unlist(intensity(be_imp)))),
                     info = paste("NA remains after", m, "imputation"))
        ## Structure unchanged
        expect_equal(length(be_imp), length(tmp))
        expect_equal(names(chromData(be_imp)), names(chromData(be)))
    }

    ## checking Gaussian specifically as it should give a warning.
    expect_warning(be_imp <- imputePeaksData(tmp, method = "gaussian"),
                   "could not fill all NAs")
    ## No more NAs
    expect_false(any(is.na(unlist(intensity(be_imp)))),
                 info = paste("NA remains after", m, "imputation"))
    ## Structure unchanged
    expect_equal(length(be_imp), length(tmp))
    expect_equal(names(chromData(be_imp)), names(chromData(be)))

    ## check loess specifically as it should give a warning.
    expect_warning(be_imp <- imputePeaksData(tmp, method = "loess"),
                   "Falling back to linear interpolation")
    ## No more NAs
    expect_false(any(is.na(unlist(intensity(be_imp)))),
                 info = paste("NA remains after", m, "imputation"))
    ## Structure unchanged
    expect_equal(length(be_imp), length(tmp))
    expect_equal(names(chromData(be_imp)), names(chromData(be)))


    ## Empty backend
    expect_silent({
        be_cd_imp <- imputePeaksData(be_empty, method = "linear")
    })
    expect_equal(length(be_cd_imp), 0L)
})

