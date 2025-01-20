test_that("peaksData works as expected", {
    peaks <- peaksData(c_full)
    backend_peaks <- peaksData(c_full@backend)
    expect_equal(peaks, backend_peaks)

    rtime_data <- peaksData(c_full, columns = "rtime", drop = TRUE)
    expect_type(rtime_data, "list")
    expect_equal(rtime_data, peaksData(c_full@backend, columns = "rtime", drop = TRUE))

    c_empty <- c_full
    c_empty@processingQueue <- list()
    expect_equal(peaksData(c_empty), peaksData(c_empty@backend))
})

test_that("peaksData replacement works as expected", {
    new_peaks <- peaksData(c_full)
    new_peaks[[1]][4, "rtime"] <- 999.9

    peaksData(c_full) <- new_peaks
    updated_peaks <- peaksData(c_full)
    expect_equal(new_peaks, updated_peaks)
})

test_that("peaksVariables works as expected", {
    vars <- peaksVariables(c_full)
    backend_vars <- peaksVariables(c_full@backend)
    expect_equal(vars, backend_vars)
    expect_true("rtime" %in% vars)
    expect_true("intensity" %in% vars)
})

test_that("rtime accessor and replacement work as expected", {
    rtime_data <- rtime(c_full)
    backend_rtime <- peaksData(c_full@backend, columns = "rtime", drop = TRUE)
    expect_equal(rtime_data, backend_rtime)

    new_rtime <- rtime(c_full)
    new_rtime[[1]][4] <- 888.8
    rtime(c_full) <- new_rtime
    expect_equal(rtime(c_full), new_rtime)
})

test_that("intensity accessor and replacement work as expected", {
    intensity_data <- intensity(c_full)
    backend_intensity <- peaksData(c_full@backend, columns = "intensity", drop = TRUE)
    expect_equal(intensity_data, backend_intensity)

    new_intensity <- intensity(c_full)
    new_intensity[[1]][1] <- 777.7
    intensity(c_full) <- new_intensity
    expect_equal(intensity(c_full), new_intensity)
})

test_that("filterPeaksData queues the correct processing step", {
    c_filtered <- filterPeaksData(c_full, variables = c("rtime"), ranges = c(12.5, 45.5))
    queue <- c_filtered@processingQueue
    expect_length(queue, 1)
    expect_equal(queue[[1]]@FUN, filterPeaksData)
    expect_equal(queue[[1]]@ARGS, list(
        variables = c("rtime"),
        ranges = c(12.5, 45.5),
        match = c("any", "all"),
        keep = TRUE
    ))

    expect_match(c_filtered@processing, "Filter: remove peaks")
})

test_that("applyProcessing applies all queued processing steps", {
    c_queued <- filterPeaksData(c_full, variables = c("rtime"),
                                ranges = c(12.5, 45.5))
    c_queued <- filterPeaksData(c_queued, variables = c("intensity"),
                                ranges = c(100, 200), keep = FALSE)
    c_applied <- applyProcessing(c_queued)
    expect_length(c_applied@processingQueue, 0)
    expect_equal(peaksData(c_applied), peaksData(c_applied@backend))
})

test_that("addProcessing adds processing steps correctly", {
    c_queued <- addProcessing(c_full, filterPeaksData,
                              variables = c("rtime"), ranges = c(12.5, 45.5))
    queue <- c_queued@processingQueue
    expect_length(queue, 1)
    expect_equal(queue[[1]]@FUN, filterPeaksData)
    expect_equal(queue[[1]]@ARGS, list(variables = c("rtime"),
                                       ranges = c(12.5, 45.5)))
    c_queued <- addProcessing(c_queued, filterPeaksData,
                              variables = c("intensity"), ranges = c(100, 200))
    queue <- c_queued@processingQueue
    expect_length(queue, 2)
})

test_that("processingChunkFactor is handled correctly in peaksData", {
    f <- processingChunkFactor(c_full)
    factorized_peaks <- peaksData(c_full, f = f)
    non_factorized_peaks <- peaksData(c_full, f = NULL)
    expect_equal(factorized_peaks, non_factorized_peaks)
})
