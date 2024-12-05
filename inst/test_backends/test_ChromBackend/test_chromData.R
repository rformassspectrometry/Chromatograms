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

test_that("backend removes NA columns", {
    chromData <- chromData(be)
    chromData$NAcol <- NA
    be2 <- backendInitialize(be, chromData)
    expect_false("NAcol" %in% names(be2@chromData))
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

test_that("dataStorage", {
    res <- dataStorage(be)
    expect_type(res, "character")
    expect_identical(length(res), length(be))
    if (!isReadOnly(be)) {
        tmp <- be
        vals <- sample(c("file", "database", "network"), length(be), replace = TRUE)
        dataStorage(tmp) <- vals
        expect_equal(dataStorage(tmp), vals)
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
})

test_that("filterChromData works", {
    res <- filterChromData(be,
                              variables = c("mz", "chromIndex"),
                              ranges = c(134, 150, 1, 1),
                              match = "any",
                              keep = FALSE)

    expect_equal(nrow(chromData(res)), 1)
    expect_equal(chromData(res)$mz, 123.3)
    expect_error(filterChromData(be,
                                 variables = c("mz"),
                                 ranges = c(134), # Wrong length
                                 match = "any",
                                 keep = TRUE),
                 "needs to be twice the length")
    expect_error(filterChromData(be,
                                 variables = c("invalid_var"),
                                 ranges = c(134, 150, 1, 1),
                                 match = "any",
                                 keep = TRUE),
                 "One or more values passed")
    res <- filterChromData(be,
                              variables = c("mz"),
                              ranges = c(134, 150),
                              match = "all",
                              keep = FALSE)

    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 123.3))

    expect_identical(filterChromData(be,
                                     match = "all",
                                     keep = TRUE),
                     be)
    expect_error(filterChromData(be,
                                 variables = c("mz"),
                                 ranges = c("a", "b"),
                                 keep = TRUE),
                 "only support filtering for numerical")
    expect_error(filterChromData(be,
                                 variables = c(200),
                                 ranges = c(134, 150),
                                 keep = TRUE),
                 "parameter needs to be of type")
    res <- filterChromData(be,
                              variables = c("mz", "chromIndex"),
                              ranges = c(134, 150, 1, 1),
                              match = "any",
                              keep = TRUE)
    expect_equal(nrow(chromData(res)), 2)
    expect_equal(chromData(res)$mz, c(112.2, 134.4))
    expect_equal(chromData(res)$chromIndex, c(1,3))
})

