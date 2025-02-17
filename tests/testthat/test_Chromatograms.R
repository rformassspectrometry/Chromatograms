
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

    expect_error(setBackend(c_mzr, backend = ChromBackendMzR()),
                 "does not support")
})

test_that("$ works correctly", {
    expect_identical(msLevel(c_full), c_full$msLevel)
    expect_identical(chromIndex(c_mzr), c_mzr$chromIndex)
    expect_identical(intensity(c_full), c_full$intensity)
    expect_identical(intensity(c_mzr), c_mzr$intensity)
    tmp <- c_full
    tmp$msLevel <- c(2L, 2L, 3L )
    expect_identical(msLevel(tmp), c(2L, 2L, 3L))
    tmp$intensity <- lapply(tmp$intensity, function(x) x + 10)
    expect_false(identical(intensity(tmp), intensity(c_full)))
})

