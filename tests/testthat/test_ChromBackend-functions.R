test_that(".valid_chrom_backend_data_storage works", {
    expect_match(.valid_chrom_backend_data_storage(c("a", NA)), "not allowed")
    expect_null(.valid_chrom_backend_data_storage(character()))
    expect_null(.valid_chrom_backend_data_storage("b"))
})

test_that(".values_match_mz works", {
    pmz <- c(12.4, 15, 3, 12.4, 3, 1234, 23, 5, 12.4, NA, 3)
    mz <- c(200, 12.4, 3)

    res <- .values_match_mz(pmz, mz)
    expect_true(all(pmz[res] %in% mz))
    expect_false(any(pmz[-res] %in% mz))

    pmz <- rev(pmz)
    res <- .values_match_mz(pmz, mz)
    expect_true(all(pmz[res] %in% mz))
    expect_false(any(pmz[-res] %in% mz))

    res <- .values_match_mz(c(NA, NA), mz)
    expect_identical(res, integer())

    res <- .values_match_mz(pmz, c(NA, 3))
    expect_true(all(pmz[res] == 3))
})

test_that("coreChromVariables works", {
    expect_equal(coreChromVariables(), .CORE_CHROM_VARIABLES)
})

test_that(".empty_peaks_data works", {
    res <- .empty_peaks_data(3, columns = c("a", "b", "c"))
    expect_true(is.list(res))
    expect_true(length(res) == 3)
    expect_equal(colnames(res[[1L]]), c("a", "b", "c"))
    expect_equal(colnames(res[[2L]]), c("a", "b", "c"))
    expect_equal(colnames(res[[3L]]), c("a", "b", "c"))
    expect_true(nrow(res[[1L]]) == 0)
    expect_true(nrow(res[[1L]]) == 0)
    expect_true(nrow(res[[1L]]) == 0)
})

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
