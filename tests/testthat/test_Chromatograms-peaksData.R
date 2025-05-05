test_that("peaksData, Chromatograms, ChrombackendMemory works as expected", {
    peaks <- peaksData(c_full)
    backend_peaks <- peaksData(c_full@backend)
    expect_equal(peaks, backend_peaks)

    rtime_data <- peaksData(c_full, columns = "rtime", drop = TRUE)
    expect_type(rtime_data, "list")
    expect_equal(rtime_data, peaksData(c_full@backend,
        columns = "rtime",
        drop = TRUE
    ))

    c_empty_queue <- c_full
    c_empty_queue@processingQueue <- list()
    expect_equal(peaksData(c_empty_queue), peaksData(c_empty_queue@backend))

    c_queued <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued@backend)))
    c_queued2 <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 1000)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued2)))

    c_queued <- filterPeaksData(c_mzr,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued@backend)))
    c_queued2 <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 1000)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued2)))
})

test_that("peaksData, Chromatogram, ChromBackendMzR works as expected", {
    peaks <- peaksData(c_mzr)
    backend_peaks <- peaksData(c_mzr@backend)
    expect_equal(peaks, backend_peaks)

    rtime_data <- peaksData(c_mzr, columns = "rtime", drop = TRUE)
    expect_type(rtime_data, "list")
    expect_equal(rtime_data, peaksData(c_mzr@backend,
        columns = "rtime",
        drop = TRUE
    ))

    c_empty_queue <- c_mzr
    c_empty_queue@processingQueue <- list()
    expect_equal(peaksData(c_empty_queue), peaksData(c_empty_queue@backend))
})

test_that("peaksData replacement works as expected", {
    new_peaks <- peaksData(c_full)
    new_peaks[[1]][4, "rtime"] <- 999.9

    peaksData(c_full) <- new_peaks
    updated_peaks <- peaksData(c_full)
    expect_equal(new_peaks, updated_peaks)

    expect_error(
        peaksData(c_mzr) <- list(data.frame()),
        "replace peaks data"
    )
})


test_that("peaksVariables works as expected", {
    vars <- peaksVariables(c_full)
    backend_vars <- peaksVariables(c_full@backend)
    expect_equal(vars, backend_vars)
    expect_true("rtime" %in% vars)
    expect_true("intensity" %in% vars)

    vars <- peaksVariables(c_mzr)
    backend_vars <- peaksVariables(c_mzr@backend)
    expect_equal(vars, backend_vars)
    expect_true("rtime" %in% vars)
    expect_true("intensity" %in% vars)
})

test_that("rtime accessor and replacement work as expected", {
    rtime_data <- rtime(c_full)
    backend_rtime <- peaksData(c_full@backend, columns = "rtime", drop = TRUE)
    expect_equal(rtime_data, backend_rtime)

    rtime_data <- rtime(c_mzr)
    backend_rtime <- peaksData(c_mzr@backend, columns = "rtime", drop = TRUE)
    expect_equal(rtime_data, backend_rtime)

    new_rtime <- rtime(c_full)
    new_rtime[[1]][4] <- 888.8
    rtime(c_full) <- new_rtime
    expect_equal(rtime(c_full), new_rtime)
    expect_error(
        rtime(c_mzr)[[1]] <- rep(1, length(rtime(c_mzr)[[1]])),
        "Cannot replace peaks data in a read-only backend"
    )
})

test_that("intensity accessor and replacement work as expected", {
    intensity_data <- intensity(c_full)
    backend_intensity <- peaksData(c_full@backend,
        columns = "intensity",
        drop = TRUE
    )
    expect_equal(intensity_data, backend_intensity)

    intensity_data <- intensity(c_mzr)
    backend_intensity <- peaksData(c_mzr@backend,
        columns = "intensity",
        drop = TRUE
    )
    expect_equal(intensity_data, backend_intensity)

    new_intensity <- intensity(c_full)
    new_intensity[[1]][1] <- 777.7
    intensity(c_full) <- new_intensity
    expect_equal(intensity(c_full), new_intensity)
    expect_error(
        intensity(c_mzr)[[1]] <- rep(1, length(intensity(c_mzr)[[1]])),
        "Cannot replace peaks data in a read-only backend"
    )
})

test_that("filterPeaksData queues the correct processing step", {
    c_filtered <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
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
