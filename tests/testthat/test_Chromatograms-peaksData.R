test_that("peaksData, Chromatograms, ChrombackendMemory works as expected", {
    peaks <- peaksData(c_full)
    backend_peaks <- peaksData(.backend(c_full))
    expect_equal(peaks, backend_peaks)

    rtime_data <- peaksData(c_full, columns = "rtime", drop = TRUE)
    expect_type(rtime_data, "list")
    expect_equal(rtime_data, peaksData(.backend(c_full),
        columns = "rtime",
        drop = TRUE
    ))

    c_empty_queue <- c_full
    c_empty_queue@processingQueue <- list()
    expect_equal(peaksData(c_empty_queue), peaksData(.backend(c_empty_queue)))

    c_queued <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    expect_false(identical(peaksData(c_queued), peaksData(.backend(c_queued))))
    c_queued2 <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 1000)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued2)))

    c_queued <- filterPeaksData(c_mzr,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    expect_false(identical(peaksData(c_queued),
                           peaksData(.backend(c_queued))))
    c_queued2 <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 1000)
    )
    expect_false(identical(peaksData(c_queued), peaksData(c_queued2)))
})

test_that("peaksData, Chromatogram, ChromBackendMzR works as expected", {
    peaks <- peaksData(c_mzr)
    backend_peaks <- peaksData(.backend(c_mzr))
    expect_equal(peaks, backend_peaks)

    rtime_data <- peaksData(c_mzr, columns = "rtime", drop = TRUE)
    expect_type(rtime_data, "list")
    expect_equal(rtime_data, peaksData(.backend(c_mzr),
        columns = "rtime",
        drop = TRUE
    ))

    c_empty_queue <- c_mzr
    c_empty_queue@processingQueue <- list()
    expect_equal(peaksData(c_empty_queue), peaksData(.backend(c_empty_queue)))
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
    backend_vars <- peaksVariables(.backend(c_full))
    expect_equal(vars, backend_vars)
    expect_true("rtime" %in% vars)
    expect_true("intensity" %in% vars)

    vars <- peaksVariables(c_mzr)
    backend_vars <- peaksVariables(.backend(c_mzr))
    expect_equal(vars, backend_vars)
    expect_true("rtime" %in% vars)
    expect_true("intensity" %in% vars)
})

test_that("rtime accessor and replacement work as expected", {
    rtime_data <- rtime(c_full)
    backend_rtime <- peaksData(.backend(c_full), columns = "rtime",
                               drop = TRUE)
    expect_equal(rtime_data, backend_rtime)

    rtime_data <- rtime(c_mzr)
    backend_rtime <- peaksData(.backend(c_mzr), columns = "rtime", drop = TRUE)
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
    backend_intensity <- peaksData(.backend(c_full),
        columns = "intensity",
        drop = TRUE
    )
    expect_equal(intensity_data, backend_intensity)

    intensity_data <- intensity(c_mzr)
    backend_intensity <- peaksData(.backend(c_mzr),
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
    queue <- .processingQueue(c_filtered)
    expect_length(queue, 1)
    expect_equal(queue[[1]]@FUN, filterPeaksData)
    expect_equal(queue[[1]]@ARGS, list(
        variables = c("rtime"),
        ranges = c(12.5, 45.5),
        match = c("any", "all"),
        keep = TRUE
    ))

    expect_match(.processing(c_filtered), "Filter: remove peaks")
})

test_that("Chromatograms, imputePeaksData works", {
    pdata <- list(
        data.frame(
            rtime = seq(12, 20, by = 0.5),
            intensity = c(123.3, 153.6, NA, 200, 210, 230, NA, 250, 260, 280,
                          300, 320, 340, 360, 380, 400, 420)
        ),
        data.frame(
            rtime = seq(45, 55, by = 0.5),  # length 21
            intensity = c(100, NA, 120, 130, 140, NA, 160, 170, 180, 190, 200,
                          210, 220, 230, 240, 250, 260, 270, 280, 290, 300)
        ),
        data.frame(
            rtime = seq(10, 18, by = 0.5),  # length 17
            intensity = c(NA, 153.6, 2354.3, 243.4, 260, NA, 280, 300, 320,
                          340, 360, 380, 400, 420, 440, 460, 480)
        )
    )

    tmp <- backendInitialize(be_empty, chromData = cdata, peaksData = pdata)
    c_tmp <- Chromatograms(tmp)

    for (meth in c("linear", "spline")) {
        c_imp <- imputePeaksData(c_tmp, method = meth, span = 0.3,
                                 sd = 1, window = 2)
        expect_s4_class(c_imp, "Chromatograms")

        expect_equal(length(peaksData(c_imp)), length(peaksData(c_tmp)))
        expect_true(all(!is.na(unlist(intensity(c_imp)))))
        expect_true(any(grepl(meth, .processing(c_imp))))
    }

    for (meth in c("gaussian", "loess")) {
        c_imp <- imputePeaksData(c_tmp, method = meth, span = 0.3,
                                 sd = 1, window = 2)
        expect_warning(c_imp <- applyProcessing(c_imp),
                       "Falling back to linear interpolation")
        expect_s4_class(c_imp, "Chromatograms")
        expect_equal(length(peaksData(c_imp)), length(peaksData(c_tmp)))
        expect_true(all(!is.na(unlist(intensity(c_imp)))))
        expect_true(any(grepl(meth, .processing(c_imp))))
    }
})
