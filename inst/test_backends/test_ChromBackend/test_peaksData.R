#' Tests for peaks variables.
test_that("peaksData", {
    res <- peaksData(be)
    expect_type(res, "list")
    expect_identical(length(res), length(be))
    expect_true(all(vapply(res, is.data.frame, logical(1))))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res), function(z) data.frame(rtime = sort(abs(rnorm(z))),
                                                            intensity = abs(rnorm(z))))
        peaksData(tmp) <- vals
        res <- peaksData(tmp)
        expect_equal(res, vals)
    }

})

test_that("peaksData() output is valid", {
    expect_true(is.null(validPeaksData(peaksData(be))))
})

test_that("peaksVariables", {
    res <- peaksVariables(be)
    expect_type(res, "character")
    expect_true(all(names(corePeaksVariables()) %in% res))

})

test_that("intensity", {
    res <- intensity(be)
    expect_true(is.list(res) || is(res, "list"))
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res), function(z) abs(rnorm(z)))
        intensity(tmp) <- vals
        res <- intensity(tmp)
        expect_equal(as.list(res), vals)
        ## test that if length is different error message
        expect_error(intensity(tmp) <- list(1:5, 1:6), "same length")
        ## test incorrect length for one element
        if (length(res) >= 2) {
            incorrect_int <- vals
            incorrect_int[[1]] <- incorrect_int[[1]][-1]
            expect_error(intensity(tmp) <- incorrect_int, "Length of ")
        }
    }
    be_empty <- ChromBackendMemory()
    expect_equal(intensity(be_empty), list())
})

test_that("rtime", {
    res <- rtime(be)
    expect_true(is.list(res) || is(res, "list"))
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- lapply(lengths(res), function(z) sort(abs(rnorm(z))))
        rtime(tmp) <- vals
        res <- rtime(tmp)
        expect_equal(as.list(res), vals)
        expect_error(rtime(tmp) <- list(1:5, 1:6), "same length")
        if (length(res) >= 2) {
            incorrect_rtime <- vals
            incorrect_rtime[[1]] <- incorrect_rtime[[1]][-1]
            expect_error(rtime(be) <- incorrect_rtime,
                         "Length of ")
        }
    }
    be_empty <- ChromBackendMemory()
    expect_equal(rtime(be_empty), list())
})


test_that("$ works", {
    res <- be$rtime
    expect_true(is.list(res) || is(res, "list"))
    expect_equal(length(res), length(be))
})

test_that("[ works", {
    be1 <- be[1]
    expect_equal(length(peaksData(be1)), 1)
    be12 <- be[1:2]
    expect_equal(length(peaksData(be12)), 2)
})



