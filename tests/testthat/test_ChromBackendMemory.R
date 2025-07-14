## Specific testing
#' - backendInitialize, ChromBackendMemory works
#' - backendMerge
#' - show
#' rest is through the generic testing in inst

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
    expect_true(all(names(.chromData(be_empty)) == names(.CORE_CHROM_VARIABLES)))
    expect_true(is.null(validChromData(.chromData(be_empty))))
    expect_true(is(.peaksData(be_empty), "list"))
    expect_true(length(.peaksData(be_empty)) == 1)
    expect_true(all(names(.peaksData(be_empty)[[1]]) == names(.CORE_PEAKS_VARIABLES)))
    expect_true(is.null(unlist(validChromData(.peaksData(be_empty)))))

    ### empty peaksData
    expect_true(is(be_cd, "ChromBackendMemory"))
    expect_true(is(.chromData(be_cd), "data.frame"))
    expect_true(nrow(.chromData(be_cd)) == nrow(cdata))
    expect_true(all(names(.chromData(be_cd)) %in% c(names(cdata), "dataOrigin")))
    expect_true(is.null(validChromData(.chromData(be_cd))))
    expect_true(all(.chromData(be_cd)$dataOrigin == NA_character_))

    expect_true(is(.peaksData(be_cd), "list"))
    expect_true(length(.peaksData(be_cd)) == nrow(.chromData(be_cd)))
    expect_true(all(vapply(.peaksData(be_cd), nrow, integer(1)) == 0))
    expect_true(all(vapply(.peaksData(be_cd), is.data.frame, logical(1))))
    expect_true(all(names(.peaksData(be_cd)[[1]]) %in% names(.CORE_PEAKS_VARIABLES)))
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
    expect_true(all(names(.peaksData(be)[[1]]) %in% names(.CORE_PEAKS_VARIABLES)))
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
    expect_true(nrow(.chromData(be_merge)) == nrow(.chromData(be_cd)) + nrow(.chromData(be)))
    expect_true(length(.peaksData(be_merge)) == length(.peaksData(be_cd)) + length(.peaksData(be)))
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
    res <- filterChromData(be,
        variables = c("mz", "chromIndex"),
        ranges = c(134, 150, 1, 1),
        match = "any",
        keep = FALSE
    )

    expect_equal(nrow(chromData(res)), 1)
    expect_equal(chromData(res)$mz, 123.3)
    expect_error(
        filterChromData(be,
            variables = c("mz"),
            ranges = c(134), # Wrong length
            match = "any",
            keep = TRUE
        ),
        "needs to be twice the length"
    )
    expect_error(
        filterChromData(be,
            variables = c("invalid_var"),
            ranges = c(134, 150, 1, 1),
            match = "any",
            keep = TRUE
        ),
        "One or more values passed"
    )
    res <- filterChromData(be,
        variables = c("mz"),
        ranges = c(134, 150),
        match = "all",
        keep = FALSE
    )

    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 123.3))

    expect_identical(
        filterChromData(be,
            match = "all",
            keep = TRUE
        ),
        be
    )
    expect_error(
        filterChromData(be,
            variables = c("mz"),
            ranges = c("a", "b"),
            keep = TRUE
        ),
        "only support filtering for numerical"
    )
    expect_error(
        filterChromData(be,
            variables = c(200),
            ranges = c(134, 150),
            keep = TRUE
        ),
        "parameter needs to be of type"
    )
    res <- filterChromData(be,
        variables = c("mz", "chromIndex"),
        ranges = c(134, 150, 1, 1),
        match = "any",
        keep = TRUE
    )
    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 134.4))
    expect_equal(chromData(res)$chromIndex, c(1, 3))

    res <- filterChromData(be,
        variables = c("mz"),
        ranges = c(999, 1000), # No matching range
        match = "any",
        keep = FALSE
    )

    expect_identical(res, be)
})


test_that("split,ChrombackendMemory works", {
    be_split <- be
    be_split$new_vars <- c("a", "b", "b")
    f <- factor(be_split$new_vars)
    res <- split(be, f)
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "ChromBackendMemory")
    expect_s4_class(res[[2L]], "ChromBackendMemory")
    expect_equal(res[[1L]]$chromIndex, c(1))
    expect_equal(res[[2L]]$chromIndex, c(2, 3))

    res <- split(be, factor(be_split$new_vars, levels = c("b", "a")))
    expect_true(is.list(res))
    expect_true(length(res) == 2)
    expect_s4_class(res[[1L]], "ChromBackendMemory")
    expect_s4_class(res[[2L]], "ChromBackendMemory")
    expect_equal(res[[2L]]$chromIndex, c(1))
    expect_equal(res[[1L]]$chromIndex, c(2, 3))
})
