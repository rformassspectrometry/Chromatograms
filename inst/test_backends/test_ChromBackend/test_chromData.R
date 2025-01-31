#' To run this unit tests from another package:
#'
#' The unit tests in this suite expect a variable `be` to be defined, which
#' has to represent an **already initialized** backend instance.

## Test chromData
test_that("chromData", {
    res <- chromData(be)
    expect_true(is.data.frame(res))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- res[sample(seq_len(nrow(res))), ]
        chromData(tmp) <- vals
        expect_equal(vals, chromData(tmp))
    }
})

test_that("chromData() output is valid", {
    expect_true(is.null(validChromData(chromData(be))))
})

test_that("chromVariables", {
    res <- chromVariables(be)
    expect_type(res, "character")
    expect_true(all(names(coreChromVariables()) %in% res))
})

test_that("chromIndex", {
    res <- chromIndex(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(1:5, length(res), replace = TRUE)
        chromIndex(tmp) <- vals
        expect_equal(chromIndex(tmp), vals)
    }
})

test_that("collisionEnergy", {
    res <- collisionEnergy(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        collisionEnergy(tmp) <- vals
        expect_equal(unname(collisionEnergy(tmp)), unname(vals))
    }
})

test_that("dataOrigin", {
    res <- dataOrigin(be)
    expect_type(res, "character")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(c("file", "database", "network"), length(be), replace = TRUE)
        dataOrigin(tmp) <- vals
        expect_equal(dataOrigin(tmp), vals)
    }
})

test_that("msLevel", {
    res <- msLevel(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(1:5, length(res), replace = TRUE)
        msLevel(tmp) <- vals
        expect_equal(msLevel(tmp), vals)
    }
})

test_that("mz", {
    res <- mz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        mz(tmp) <- vals
        expect_equal(unname(mz(tmp)), unname(vals))
    }
})

test_that("mzMax", {
    res <- mzMax(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        mzMax(tmp) <- vals
        expect_equal(unname(mzMax(tmp)), unname(vals))
    }
})

test_that("mzMin", {
    res <- mzMin(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        mzMin(tmp) <- vals
        expect_equal(unname(mzMin(tmp)), unname(vals))
    }
})

test_that("precursorMz", {
    res <- precursorMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        precursorMz(tmp) <- vals
        expect_equal(unname(precursorMz(tmp)), unname(vals))
    }
})

test_that("precursorMzMin", {
    res <- precursorMzMin(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        precursorMzMin(tmp) <- vals
        expect_equal(unname(precursorMzMin(tmp)), unname(vals))
    }
})

test_that("precursorMzMax", {
    res <- precursorMzMax(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        precursorMzMax(tmp) <- vals
        expect_equal(unname(precursorMzMax(tmp)), unname(vals))
    }
})

test_that("productMz", {
    res <- productMz(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        productMz(tmp) <- vals
        expect_equal(unname(productMz(tmp)), unname(vals))
    }
})

test_that("productMzMin", {
    res <- productMzMin(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        productMzMin(tmp) <- vals
        expect_equal(unname(productMzMin(tmp)), unname(vals))
    }
})

test_that("productMzMax", {
    res <- productMzMax(be)
    expect_type(res, "double")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- abs(rnorm(length(be)))
        productMzMax(tmp) <- vals
        expect_equal(unname(productMzMax(tmp)), unname(vals))
    }
})

## Other
test_that("isEmpty", {
    res <- isEmpty(be)
    expect_type(res, "logical")
    expect_identical(length(res), length(be))
    expect_identical(lengths(be) == 0, res)
})

test_that("lengths", {
    res <- lengths(be)
    expect_type(res, "integer")
    expect_identical(length(res), length(be))
    expect_identical(res, lengths(rtime(be)))
})

test_that("$ works", {
    res <- be$msLevel
    expect_true(is.integer(res))
    expect_equal(length(res), length(be))

    ## Error if spectra variable not available
    expect_error(be$doesnt_exist)
})

test_that("[ works", {
    be1 <- be[1]
    expect_equal(nrow(chromData(be1)), 1)
    expect_equal(chromData(be1)$chromIndex, 1)
    be12 <- be[1:2]
    expect_equal(nrow(chromData(be12)), 2)
    expect_equal(chromData(be12)$chromIndex, c(1, 2))

    # Handles empty integer
    be0 <- be[integer()]
    expect_equal(nrow(chromData(be0)), 0)
    expect_equal(be0, be_empty)
})

test_that("[[ works", {
    expect_error(be[["doesnt_exist"]], "The requested variable")
    expect_error(be[["msLevel", "no"]], "is not supported")
    expect_equal(be[["msLevel"]], msLevel(be))
    if (!isReadOnly(be)) {
        be2 <- be
        be2[["msLevel"]] <- c(1L,2L,3L)
        expect_false(all(be2[["msLevel"]] == msLevel(be)))
    }
})

