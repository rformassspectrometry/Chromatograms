test_that("fillCoreChromVariables works", {
    x <- data.frame(a = numeric(), b = character())
    res <- fillCoreChromVariables(x)
    expect_true(is.data.frame(res))
    expect_true(all(c("a", "b") %in% colnames(res)))
    expect_true(all(names(.CORE_CHROM_VARIABLES) %in%
                        c(colnames(res), "intensity", "rtime")))
    expect_true(nrow(res) == 0)

    x <- data.frame(msLevel = c(1L, 2L, 1L), other_col = "a")
    res <- fillCoreChromVariables(x)
    expect_true(is.data.frame(res))
    expect_true(all(c("msLevel", "other_col") %in% colnames(res)))
    expect_true(all(names(.CORE_CHROM_VARIABLES) %in%
                        c(colnames(res), "intensity", "rtime")))
    expect_true(nrow(res) == 3L)

    cv <- .CORE_CHROM_VARIABLES[!names(.CORE_CHROM_VARIABLES) %in%
                                    c("intensity", "rtime", "msLevel")]
    for (i in seq_along(cv)) {
        expect_true(is(res[, names(cv)[i]], cv[i]))
        expect_true(all(is.na(res[, names(cv)[i]])))
    }
    expect_equal(res$msLevel, c(1L, 2L, 1L))

    x <- fillCoreChromVariables(x)
    expect_equal(x, fillCoreChromVariables(x))
})

test_that("core variables functions works", {
    expect_equal(coreChromVariables(), .CORE_CHROM_VARIABLES)
    expect_equal(corePeaksVariables(), .CORE_PEAKS_VARIABLES)
})

test_that("validChromData works", {
    x <- data.frame()
    res <- validChromData(x)
    expect_true(length(res) == 0)

    x <- data.frame(msLevel = c(1L, 2L), other_col = "a")
    res <- validChromData(x)
    expect_true(length(res) == 0)

    x$mz <- "a"
    x$collisionEnergy <- "30ev"
    expect_error(validChromData(x, error = TRUE), "wrong data type")
    res <- validChromData(x, error = FALSE)
    expect_true(is.character(res))
    expect_true(length(res) == 2L)
})

test_that("validPeaksData works", {
    x <- matrix()
    expect_error(validPeaksData(x), "list")
    x <- list(data.frame(rtime = numeric(), intensity = numeric()))
    res <- validPeaksData(x)
    expect_true(is.null(res))
    #checks order
    x <- list(data.frame(intensity = numeric(), rtime = numeric()))
    expect_error(validPeaksData(x), "Columns should be in the order")

    x <- list(
        data.frame(rtime = c(10.0, 12.0), intensity = c(200, 150)),
        data.frame(rtime = c(30.1, 31.2), intensity = c(110, 90), other_col= c("test", "test"))
    )
    expect_error(validPeaksData(x),
                 "same order")
})

