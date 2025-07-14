test_that("addProcessing adds processing steps correctly", {
    c_queued <- addProcessing(c_full, filterPeaksData,
        variables = c("rtime"), ranges = c(12.5, 45.5)
    )
    queue <- .processingQueue(c_queued)
    expect_length(queue, 1)
    expect_equal(queue[[1]]@FUN, filterPeaksData)
    expect_equal(queue[[1]]@ARGS, list(
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    ))
    c_queued <- addProcessing(c_queued, filterPeaksData,
        variables = c("intensity"), ranges = c(100, 200)
    )
    queue <- .processingQueue(c_queued)
    expect_length(queue, 2)

    ## Missing FUNs
    c_queued <- addProcessing(c_mzr,
        variables = c("rtime"), ranges = c(12.5, 45.5)
    )
    expect_equal(c_queued, c_mzr)
})

test_that("processingChunkFactor is handled correctly in peaksData", {
    f <- processingChunkFactor(c_full, chunkSize = 2)
    factorized_peaks <- peaksData(c_full, f = f)
    non_factorized_peaks <- peaksData(c_full, f = f())
    expect_equal(factorized_peaks, non_factorized_peaks)

    f <- processingChunkFactor(c_mzr, chunkSize = 100)
    factorized_peaks <- peaksData(c_mzr, f = f)
    non_factorized_peaks <- peaksData(c_mzr, f = f())
    expect_equal(factorized_peaks, non_factorized_peaks)
})

test_that("applyProcessing applies all queued processing steps", {
    c_queued <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    c_queued <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 200), keep = FALSE
    )
    c_applied <- applyProcessing(c_queued)
    expect_length(.processingQueue(c_applied), 0)
    expect_equal(peaksData(c_applied), peaksData(.backend(c_applied)))

    ## empty queue
    c_queued@processingQueue <- list()
    c_applied <- applyProcessing(c_queued)
    expect_equal(c_applied, c_queued)

    c_queued <- filterPeaksData(c_mzr,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    expect_error(applyProcessing(c_queued), "Cannot apply")
})
