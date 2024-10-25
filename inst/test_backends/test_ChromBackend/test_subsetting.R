#' split:
test_that("split", {
    res <- split(be, f = seq_along(be))
    expect_true(is.list(res))
    expect_equal(length(res), length(be))
    for (i in seq_along(be)) {
        expect_s4_class(res[[i]], class(be)[1L])
        expect_true(is.null(validChromData(chromData(res[[i]]))))
        expect_true(is.null(validPeaksData(peaksData(res[[i]]))))
        expect_true(length(res[[i]]) == 1L)
    }
})

#'[ , any order, duplication.
test_that("[", {
    set.seed(123)
    ## random order
    idx <- sample(seq_along(be))
    res <- be[idx]
    expect_equal(length(res), length(be))
    expect_true(is.null(validChromData(chromData(res))))
    expect_true(is.null(validPeaksData(peaksData(res))))
    for (i in seq_along(idx)) {
        a <- chromData(res[i])
        b <- chromData(be[idx[i]])
        expect_equal(a, b)
    }

    ## duplication
    res <- be[c(1, 1, 1)]
    expect_equal(length(res), 3L)
    expect_true(is.null(validChromData(chromData(res))))
    expect_true(is.null(validPeaksData(peaksData(res))))
    a <- chromData(be[1L])
    rownames(a) <- NULL
    b <- chromData(res[1L])
    rownames(b) <- NULL
    expect_equal(a, b)
    b <- chromData(res[2L])
    rownames(b) <- NULL
    expect_equal(a, b)
    b <- chromData(res[3L])
    rownames(b) <- NULL
    expect_equal(a, b)

    ## Out of range should throw error
    expect_error(be[0], "index")

    ## integer(0) should return an empty object
    res <- be[integer()]
    expect_s4_class(res, class(be)[1L])
    expect_true(length(res) == 0L)
})

#' test if any eventually implemented method yields the same result as the
#' default implementation
test_that("filterDataOrigin", {
    ref <- getMethod("filterDataOrigin", "ChromBackend")
    org <- unique(dataOrigin(be))[1L]
    if (!is.na(org)) {
        a <- ref(be, org)
        b <- filterDataOrigin(be, org)
        a <- chromData(a)
        b <- chromData(b)
        expect_equal(a, b)
    }
})

test_that("filterMsLevel", {
    expect_equal(be, filterMsLevel(be))
    ref <- getMethod("filterMsLevel", "ChromBackend")
    org <- unique(msLevel(be))[1L]
    if (!is.na(org)) {
        a <- ref(be, org)
        b <- filterMsLevel(be, org)
        a <- chromData(a)
        b <- chromData(b)
        expect_equal(a, b)
    }
})

test_that("filterMzRange", {
    expect_equal(be, filterMzRange(be))
    ref <- getMethod("filterMzRange", "ChromBackend")
    mz <- range(mz(be))
    if (!is.na(mz[1L])) {
        a <- ref(be, mz)
        b <- filterMzRange(be, mz)
        a <- chromData(a)
        b <- chromData(b)
        expect_equal(a, b)
    }
})

test_that("filterMzValues", {
    expect_equal(be, filterMzValues(be))
    ref <- getMethod("filterMzValues", "ChromBackend")
    mz <- mz(be)
    if (!is.na(mz[1L])) {
        a <- ref(be, mz)
        b <- filterMzValues(be, mz)
        a <- chromData(a)
        b <- chromData(b)
        expect_equal(a, b)
    }
})

