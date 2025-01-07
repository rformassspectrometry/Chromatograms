

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
})

test_that("show, Chromatograms works", {
    expect_output(show(c_full), "ChromBackendMemory")
    res <- c_full
    res@processing <- c("a", "b", "c", "d")
    expect_output(show(res), "1 more processings")
    res@processingQueue <- list("a", "b", "c", "d")
    expect_output(show(res), "4 processing step")
    show(res)

})

