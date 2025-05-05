test_that("Chromatograms works", {
    ## empty object
    expect_true(is(c_empty@backend, "ChromBackendMemory"))
    expect_true(is(c_empty, "Chromatograms"))
    expect_true(c_empty@processingChunkSize == Inf)
    expect_true(c_empty@version == "0.1")
    expect_identical(c_empty@processingQueue, list())

    ## object with backend
    expect_true(is(c_full@backend, "ChromBackendMemory"))
    expect_true(is(c_full, "Chromatograms"))
    expect_true(c_full@processingChunkSize == Inf)
    expect_true(c_full@version == "0.1")
    expect_identical(c_full@processingQueue, list())

    expect_equal(processingChunkSize(c_full), Inf)
    c_chunk <- c_full
    processingChunkSize(c_chunk) <- 2
    expect_equal(processingChunkSize(c_chunk), 2)
    expect_equal(levels(processingChunkFactor(c_chunk)), c("1", "2"))

    ## method with Spectra works
    c_sp <- Chromatograms(s[1:2])
    expect_true(is(c_sp@backend, "ChromBackendSpectra"))
    expect_true(is(c_sp, "Chromatograms"))
    expect_true(c_sp@processingChunkSize == Inf)
    expect_true(c_sp@version == "0.1")
    expect_identical(c_sp@processingQueue, list())
})

test_that("show, Chromatograms - ChromBackendMemory works", {
    expect_output(show(c_full), "ChromBackendMemory")
    res <- c_full
    res@processing <- c("a", "b", "c", "d")
    expect_output(show(res), "1 more processings")
    res@processingQueue <- list("a", "b", "c", "d")
    expect_output(show(res), "4 processing step")
})

test_that("show, Chromatograms - ChromBackendMzR works", {
    expect_output(show(c_mzr), "ChromBackendMzR")
    res <- c_mzr
    res@processing <- c("a", "b", "c", "d")
    expect_output(show(res), "1 more processings")
    res@processingQueue <- list("a", "b", "c", "d")
    expect_output(show(res), "4 processing step")
})


test_that("setBackend works correctly", {
    c_mzr_new <- setBackend(c_mzr, backend = ChromBackendMemory())
    expect_s4_class(c_mzr_new@backend, "ChromBackendMemory")
    expect_identical(chromData(c_mzr_new), chromData(c_mzr))
    expect_identical(peaksData(c_mzr_new), peaksData(c_mzr))
    expect_identical(c_mzr_new@backend@peaksData, peaksData(c_mzr))

    processingChunkSize(c_mzr) <- 100
    f <- processingChunkFactor(c_mzr)
    expect_true(length(levels(f)) > 1)
    c_mzr_new <- setBackend(c_mzr, backend = ChromBackendMemory(), f = f)
    expect_s4_class(c_mzr_new@backend, "ChromBackendMemory")
    expect_identical(chromData(c_mzr_new), chromData(c_mzr))
    expect_identical(peaksData(c_mzr_new), peaksData(c_mzr))
    expect_identical(c_mzr_new@backend@peaksData, peaksData(c_mzr))

    expect_error(
        setBackend(c_mzr, backend = ChromBackendMzR()),
        "does not support"
    )
})

test_that("$ works correctly", {
    expect_identical(msLevel(c_full), c_full$msLevel)
    expect_identical(chromIndex(c_mzr), c_mzr$chromIndex)
    expect_identical(intensity(c_full), c_full$intensity)
    expect_identical(intensity(c_mzr), c_mzr$intensity)
    tmp <- c_full
    tmp$msLevel <- c(2L, 2L, 3L)
    expect_identical(msLevel(tmp), c(2L, 2L, 3L))
    tmp$intensity <- lapply(tmp$intensity, function(x) x + 10)
    expect_false(identical(intensity(tmp), intensity(c_full)))
})

test_that("[ works correctly", {
    c_sub <- c_full[1:2]
    expect_true(is(c_sub, "Chromatograms"))
    expect_equal(nrow(chromData(c_sub)), 2)
    expect_equal(length(peaksData(c_sub)), 2)
    expect_error(c_full[1:2, 1], "by columns is not")

    c_sub <- c_full[1]
    expect_equal(c_sub, c_sub[])
})

test_that("[[ works properly", {
    expect_error(c_full[[1]], "character")
    expect_error(c_full[["test", 1]], "not supported")
    expect_error(c_full[["test"]], "No variable")
    expect_equal(c_full[["msLevel"]], msLevel(c_full))
    expect_equal(c_full[["msLevel"]], c_full@backend[["msLevel"]])

    ## replace
    expect_error(c_full[[1]] <- 1, "character defining the chromatogram")
    expect_error(c_full[["msLevel", 4]] <- 1, "not supported")
    expect_error(c_full[["test"]] <- 1, "No variable")
    repet <- c_full
    repet[["msLevel"]] <- rep(2L, length(repet))
    expect_false(identical(msLevel(repet), msLevel(c_full)))
})

test_that("factorize() works", {
    tmp <- c_sp
    tmp$msLevel <- c(1L, 2L, 3L)
    idx_before <- chromSpectraIndex(tmp@backend)
    tmp <- factorize(tmp)
    idx_after <- chromSpectraIndex(tmp@backend)
    expect_false(identical(idx_before, idx_after))
})
