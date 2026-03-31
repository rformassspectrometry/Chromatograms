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
                                 sd = 1, window = 2, extrapolate = TRUE)
        expect_s4_class(c_imp, "Chromatograms")

        expect_equal(length(peaksData(c_imp)), length(peaksData(c_tmp)))
        expect_true(all(!is.na(unlist(intensity(c_imp)))))
        expect_true(any(grepl(meth, .processing(c_imp))))
    }

    for (meth in c("gaussian", "loess")) {
        c_imp <- imputePeaksData(c_tmp, method = meth, span = 0.3,
                                 sd = 1, window = 2, extrapolate = TRUE)
        expect_warning(c_imp <- applyProcessing(c_imp,
                                                BPPARAM = SerialParam()),
                       "Falling back to linear interpolation")
        expect_s4_class(c_imp, "Chromatograms")
        expect_equal(length(peaksData(c_imp)), length(peaksData(c_tmp)))
        expect_true(all(!is.na(unlist(intensity(c_imp)))))
        expect_true(any(grepl(meth, .processing(c_imp))))
    }

    ## Test extrapolate = FALSE (default)
    ## Interior NAs should be interpolated, leading/trailing NAs remain
    pdata_edge <- list(
        data.frame(
            rtime = 1:10,
            intensity = c(NA, NA, 100, NA, 140, 160, 180, NA, NA, NA)
        )
    )
    tmp_edge <- backendInitialize(be_empty, chromData = cdata[1, , drop = FALSE],
                                  peaksData = pdata_edge)
    c_edge <- Chromatograms(tmp_edge)

    ## extrapolate = FALSE (default) - leading/trailing NAs remain, interior filled
    c_imp_interp <- imputePeaksData(c_edge, method = "linear")
    ints <- intensity(c_imp_interp)[[1]]
    expect_true(is.na(ints[1]))
    expect_true(is.na(ints[2]))
    expect_true(is.na(ints[8]))
    expect_true(is.na(ints[9]))
    expect_true(is.na(ints[10]))
    # interior NA at position 4 should be interpolated
    expect_false(is.na(ints[4]))
    # all interior values (positions 3-7) should have no NAs
    expect_false(any(is.na(ints[3:7])))

    ## extrapolate = TRUE - no NAs remain
    c_imp_extrap <- imputePeaksData(c_edge, method = "linear", extrapolate = TRUE)
    expect_false(any(is.na(unlist(intensity(c_imp_extrap)))))

    ## Check processing log mentions extrapolation when enabled
    expect_true(grepl("with extrapolation", .processing(c_imp_extrap)))
})

test_that("peakBoundary returns per-chromatogram matrix.", {
    cdata <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100.0, 200.0),
        dataOrigin = c("mem1", "mem1")
    )
    pdata <- list(
        data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
                   intensity = c(100, 250, 400, 300, 150)),
        data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
                   intensity = c(0, 10, 50, 100, 50, 10, 0))
    )
    chr <- Chromatograms(ChromBackendMemory(), chromData = cdata,
                         peaksData = pdata)
    tmp <- peakBoundary(chr)
    expect_true(is.matrix(tmp))
    expect_equal(nrow(tmp), length(chr))
    expect_equal(colnames(tmp), c("left_boundary", "right_boundary"))
    expect_true(all(tmp[, "left_boundary"] <= tmp[, "right_boundary"],
                    na.rm = TRUE))
})

test_that("peakBoundary returns correct boundaries for clean symmetric peak.", {
    cdata_pb <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_pb <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(0, 10, 50, 100, 50, 10, 0)
    ))
    chr_pb <- Chromatograms(ChromBackendMemory(), chromData = cdata_pb,
                            peaksData = pdata_pb)
    tmp <- peakBoundary(chr_pb)
    expect_equal(unname(tmp[1, "left_boundary"]), 1)
    expect_equal(unname(tmp[1, "right_boundary"]), 7)
})

test_that("peakBoundary handles empty chromatogram.", {
    cdata_e <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_e <- list(data.frame(rtime = numeric(0), intensity = numeric(0)))
    chr_empty <- Chromatograms(ChromBackendMemory(), chromData = cdata_e,
                               peaksData = pdata_e)
    tmp <- peakBoundary(chr_empty)
    expect_true(all(is.na(tmp)))

    ## Also n < 3 returns NA
    pdata_2 <- list(data.frame(rtime = c(1, 2), intensity = c(10, 20)))
    chr_2 <- Chromatograms(ChromBackendMemory(), chromData = cdata_e,
                           peaksData = pdata_2)
    tmp2 <- peakBoundary(chr_2)
    expect_true(all(is.na(tmp2)))
})

test_that("peakBoundary handles all-NA intensities.", {
    cdata_na <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_na <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5),
        intensity = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
    ))
    chr_allna <- Chromatograms(ChromBackendMemory(), chromData = cdata_na,
                               peaksData = pdata_na)
    tmp <- peakBoundary(chr_allna)
    expect_true(all(is.na(tmp)))
})

test_that("peakBoundary handles all-zero intensities.", {
    cdata_z <- data.frame(msLevel = c(1L, 1L), mz = c(100, 200),
                          dataOrigin = c("s1", "s1"))
    pdata_z <- list(
        data.frame(rtime = c(1, 2, 3, 4, 5),
                   intensity = c(0, 0, 0, 0, 0)),
        data.frame(rtime = c(1, 2, 3, 4, 5, 6, 7),
                   intensity = c(0, 0, 0, 0, 0, 0, 0))
    )
    chr_zero <- Chromatograms(ChromBackendMemory(), chromData = cdata_z,
                              peaksData = pdata_z)
    tmp <- peakBoundary(chr_zero)
    expect_true(all(is.na(tmp)))
})

test_that("peakBoundary threshold parameter works.", {
    cdata_t <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_t <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7),
        intensity = c(5, 20, 60, 100, 60, 20, 5)
    ))
    chr_t <- Chromatograms(ChromBackendMemory(), chromData = cdata_t,
                           peaksData = pdata_t)
    tmp <- peakBoundary(chr_t, threshold = 0.1)
    expect_true(is.matrix(tmp))
    expect_false(any(is.na(tmp)))
    tmp_05 <- peakBoundary(chr_t, threshold = 0.05)
    expect_true(is.matrix(tmp_05))
    expect_false(any(is.na(tmp_05)))
})

test_that("peakBoundary correctly isolates tallest peak in multi-peak data.", {
    cdata_m <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_m <- list(data.frame(
        rtime = 1:9,
        intensity = c(0, 5, 10, 5, 0, 3, 20, 3, 0)
    ))
    chr_m <- Chromatograms(ChromBackendMemory(), chromData = cdata_m,
                           peaksData = pdata_m)
    tmp <- peakBoundary(chr_m)
    ## Should return boundaries for the tallest peak (rtime 5-9 region)
    expect_equal(unname(tmp[1, "left_boundary"]), 5)
    expect_equal(unname(tmp[1, "right_boundary"]), 9)
})

test_that("peakBoundary handles NA-adjacent boundaries.", {
    cdata_na_adj <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_na_adj <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        intensity = c(NA, NA, 50, 100, 500, 100, 50, NA, NA, 5)
    ))
    chr_na_adj <- Chromatograms(ChromBackendMemory(), chromData = cdata_na_adj,
                                peaksData = pdata_na_adj)
    expect_no_error(peakBoundary(chr_na_adj))
    tmp <- peakBoundary(chr_na_adj)
    expect_true(is.matrix(tmp))
    expect_false(any(is.na(tmp)))
})

test_that("peakBoundary validates threshold parameter.", {
    cdata_v <- data.frame(msLevel = 1L, mz = 100, dataOrigin = "s1")
    pdata_v <- list(data.frame(
        rtime = c(1, 2, 3, 4, 5),
        intensity = c(10, 50, 100, 50, 10)
    ))
    chr_v <- Chromatograms(ChromBackendMemory(), chromData = cdata_v,
                           peaksData = pdata_v)
    expect_error(peakBoundary(chr_v, threshold = "a"),
                 "single non-missing numeric")
    expect_error(peakBoundary(chr_v, threshold = c(0.1, 0.2)),
                 "single non-missing numeric")
    expect_error(peakBoundary(chr_v, threshold = NA),
                 "single non-missing numeric")
    expect_error(peakBoundary(chr_v, threshold = -0.1),
                 "must be >= 0 and < 1")
    expect_error(peakBoundary(chr_v, threshold = 1),
                 "must be >= 0 and < 1")
    expect_error(peakBoundary(chr_v, threshold = 2),
                 "must be >= 0 and < 1")
    ## baselineThreshold validation
    expect_error(peakBoundary(chr_v, baselineThreshold = "a"),
                 "baselineThreshold.*single non-missing numeric")
    expect_error(peakBoundary(chr_v, baselineThreshold = -0.1),
                 "baselineThreshold.*must be >= 0 and < 1")
    ## baselineQuantile validation
    expect_error(peakBoundary(chr_v, baselineQuantile = "a"),
                 "baselineQuantile.*single non-missing numeric")
    expect_error(peakBoundary(chr_v, baselineQuantile = 1.5),
                 "baselineQuantile.*must be >= 0 and <= 1")
})

## ---- compareChromatograms tests ----

test_that("compareChromatograms default returns correct array structure.", {
    res <- compareChromatograms(c_full)
    expect_true(is.array(res))
    expect_equal(dim(res), c(3L, 3L, 2L))
    expect_equal(dimnames(res)[[3L]], c("score", "n_peaks"))
    ## Scores layer: symmetric, diagonal = 1
    expect_equal(diag(res[, , 1L]), c(1, 1, 1))
    expect_equal(res[1, 2, 1L], res[2, 1, 1L])
    expect_equal(res[1, 3, 1L], res[3, 1, 1L])
    expect_equal(res[2, 3, 1L], res[3, 2, 1L])
    expect_equal(unname(res[1, 3, 1L]), 1)
    expect_true(is.na(res[1, 2, 1L]))
    expect_true(is.na(res[2, 3, 1L]))
    ## Count layer: diagonal > 0, non-overlapping = 0
    expect_true(all(diag(res[, , 2L]) > 0))
    expect_equal(unname(res[1, 2, 2L]), 0)
    expect_equal(unname(res[2, 3, 2L]), 0)
    expect_null(rownames(res))
    expect_null(colnames(res))
})

test_that("compareChromatograms default, empty object returns 0x0x2 array.", {
    res <- compareChromatograms(c_empty)
    expect_true(is.array(res))
    expect_equal(dim(res), c(0L, 0L, 2L))
})

test_that("compareChromatograms default, single chromatogram returns 1x1x2.", {
    res <- compareChromatograms(c_full[1])
    expect_true(is.array(res))
    expect_equal(dim(res), c(1L, 1L, 2L))
    expect_equal(unname(res[1, 1, 1L]), 1)
    expect_true(res[1, 1, 2L] > 0)
})

test_that("compareChromatograms method = 'spearman' and 'kendall' work.", {
    res_s <- compareChromatograms(c_full, method = "spearman")
    expect_true(is.array(res_s))
    expect_equal(dim(res_s), c(3L, 3L, 2L))
    expect_equal(diag(res_s[, , 1L]), c(1, 1, 1))
    ## Identical chromatograms → 1 regardless of method
    expect_equal(unname(res_s[1, 3, 1L]), 1)
    ## Non-overlapping → NA
    expect_true(is.na(res_s[1, 2, 1L]))

    res_k <- compareChromatograms(c_full, method = "kendall")
    expect_true(is.array(res_k))
    expect_equal(diag(res_k[, , 1L]), c(1, 1, 1))
    expect_equal(unname(res_k[1, 3, 1L]), 1)
})

test_that("compareChromatograms with custom FUN produces expected values.", {
    cosine_sim <- function(x, y, ...) {
        sum(x * y, na.rm = TRUE) /
            (sqrt(sum(x^2, na.rm = TRUE)) * sqrt(sum(y^2, na.rm = TRUE)))
    }
    res <- compareChromatograms(c_full, FUN = cosine_sim)
    expect_true(is.array(res))
    expect_equal(dim(res), c(3L, 3L, 2L))
    expect_equal(diag(res[, , 1L]), c(1, 1, 1))
    expect_equal(unname(res[1, 3, 1L]), 1)
})

test_that("compareChromatograms returns NA for non-overlapping chromatograms.", {
    cdata_no <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100, 200),
        dataOrigin = c("s1", "s2")
    )
    pdata_no <- list(
        data.frame(rtime = c(1, 2, 3, 4), intensity = c(10, 50, 100, 50)),
        data.frame(rtime = c(10, 11, 12, 13), intensity = c(20, 80, 200, 80))
    )
    chr_no <- Chromatograms(ChromBackendMemory(), chromData = cdata_no,
                            peaksData = pdata_no)
    res <- compareChromatograms(chr_no)
    expect_true(is.na(res[1, 2, 1L]))
    expect_true(is.na(res[2, 1, 1L]))
    expect_equal(diag(res[, , 1L]), c(1, 1))
})

test_that("compareChromatograms returns NA when overlap < minPeaks.", {
    cdata_short <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100, 200),
        dataOrigin = c("s1", "s2")
    )
    pdata_short <- list(
        data.frame(rtime = 1, intensity = 50),
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 50, 100))
    )
    chr_short <- Chromatograms(ChromBackendMemory(), chromData = cdata_short,
                               peaksData = pdata_short)
    res <- compareChromatograms(chr_short)
    expect_true(is.na(res[1, 2, 1L]))
    expect_equal(unname(res[1, 1, 1L]), 1)
    expect_equal(unname(res[2, 2, 1L]), 1)
})

test_that("compareChromatograms computes correct Pearson value for known data.", {
    ## Two chromatograms with overlapping RT and known correlation
    cdata_known <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100, 200),
        dataOrigin = c("s1", "s2")
    )
    pdata_known <- list(
        data.frame(rtime = c(1, 2, 3, 4, 5),
                   intensity = c(10, 20, 30, 20, 10)),
        data.frame(rtime = c(1, 2, 3, 4, 5),
                   intensity = c(5, 15, 25, 15, 5))
    )
    chr_known <- Chromatograms(ChromBackendMemory(), chromData = cdata_known,
                               peaksData = pdata_known)
    res <- compareChromatograms(chr_known)
    ## Same RT grid → no interpolation needed; these are perfectly correlated
    expect_equal(unname(res[1, 2, 1L]), 1, tolerance = 1e-10)
    expect_equal(unname(res[1, 2, 2L]), 5)  ## 5 overlapping RT points
})

test_that("compareChromatograms minPeaks filters low-overlap pairs.", {
    ## chr_known has 5 overlapping RT points
    cdata_known <- data.frame(
        msLevel = c(1L, 1L),
        mz = c(100, 200),
        dataOrigin = c("s1", "s2")
    )
    pdata_known <- list(
        data.frame(rtime = c(1, 2, 3, 4, 5),
                   intensity = c(10, 20, 30, 20, 10)),
        data.frame(rtime = c(1, 2, 3, 4, 5),
                   intensity = c(5, 15, 25, 15, 5))
    )
    chr_known <- Chromatograms(ChromBackendMemory(), chromData = cdata_known,
                               peaksData = pdata_known)
    ## Default minPeaks = 4: 5 points ≥ 4 → score computed
    res <- compareChromatograms(chr_known)
    expect_false(is.na(res[1, 2, 1L]))
    ## minPeaks = 6: 5 points < 6 → score is NA, but count is still 5
    res_strict <- compareChromatograms(chr_known, minPeaks = 6L)
    expect_true(is.na(res_strict[1, 2, 1L]))
    expect_equal(unname(res_strict[1, 2, 2L]), 5)
    ## minPeaks = 2: always compute when ≥ 2 points overlap
    res_loose <- compareChromatograms(chr_known, minPeaks = 2L)
    expect_false(is.na(res_loose[1, 2, 1L]))
})

test_that("compareChromatograms labels parameter sets row/col names.", {
    ## 'mz' column has unique values in c_full
    res <- compareChromatograms(c_full, labels = "mz")
    expect_equal(dimnames(res)[[1L]], as.character(chromData(c_full)[["mz"]]))
    expect_equal(dimnames(res)[[2L]], as.character(chromData(c_full)[["mz"]]))
    ## NULL labels → no names on first two dims
    res_no <- compareChromatograms(c_full, labels = NULL)
    expect_null(dimnames(res_no)[[1L]])
    expect_null(dimnames(res_no)[[2L]])
    ## Non-existent column → error
    expect_error(compareChromatograms(c_full, labels = "nonexistent"),
                 "not found")
    ## Non-unique column → error
    c_dup <- c_full
    c_dup$msLevel <- rep(1L, length(c_dup))
    expect_error(compareChromatograms(c_dup, labels = "msLevel"),
                 "duplicated")
    ## Non-character labels → error
    expect_error(compareChromatograms(c_full, labels = 42),
                 "single character string")
})

test_that("compareChromatograms works with ChromBackendMzR.", {
    res <- compareChromatograms(c_mzr)
    expect_true(is.array(res))
    expect_equal(dim(res)[1L], length(c_mzr))
    expect_equal(dim(res)[2L], length(c_mzr))
    expect_equal(dim(res)[3L], 2L)
    expect_equal(diag(res[, , 1L]), rep(1, length(c_mzr)))
    expect_true(is.numeric(res))
})

test_that("compareChromatograms per-group via split() works.", {
    ## The 'by' parameter is no longer supported; users split beforehand.
    c_grp <- c_full
    c_grp$group <- c("A", "B", "A")
    grp_list <- split(c_grp, c_grp$group)
    res <- lapply(grp_list, compareChromatograms)
    expect_equal(sort(names(res)), c("A", "B"))
    ## Group "A" has 2 identical chromatograms → off-diagonal score = 1
    expect_equal(dim(res[["A"]]), c(2L, 2L, 2L))
    expect_equal(unname(res[["A"]][1, 2, 1L]), 1)
    ## Group "B" has 1 chromatogram → 1x1x2 array
    expect_equal(dim(res[["B"]]), c(1L, 1L, 2L))
    expect_equal(unname(res[["B"]][1, 1, 1L]), 1)
})

test_that("compareChromatograms minPeaks works for x,y comparison.", {
    ## chr_known has 5 overlapping RT points (defined above in minPeaks test)
    cdata_known <- data.frame(
        msLevel = c(1L, 1L), mz = c(100, 200), dataOrigin = c("s1", "s2")
    )
    pdata_known <- list(
        data.frame(rtime = c(1, 2, 3, 4, 5), intensity = c(10, 20, 30, 20, 10)),
        data.frame(rtime = c(1, 2, 3, 4, 5), intensity = c(5, 15, 25, 15, 5))
    )
    chr_known <- Chromatograms(ChromBackendMemory(), chromData = cdata_known,
                               peaksData = pdata_known)
    ## strict threshold blocks score but count is still recorded
    res <- compareChromatograms(chr_known[1], chr_known[2], minPeaks = 6L)
    expect_true(is.na(res[1, 1, 1L]))
    expect_equal(unname(res[1, 1, 2L]), 5)
})

## ---- compareChromatograms x, y tests ----

test_that("compareChromatograms(x, y) returns n x m x 2 array.", {
    ## c_full has 3 chromatograms; compare first 2 vs last 1
    res <- compareChromatograms(c_full[1:2], c_full[3])
    expect_true(is.array(res))
    expect_equal(dim(res), c(2L, 1L, 2L))
    ## Chromatograms 1 and 3 are identical
    expect_equal(unname(res[1, 1, 1L]), 1)
    ## Chromatograms 2 and 3 have non-overlapping RT
    expect_true(is.na(res[2, 1, 1L]))
    ## Overlap count for identical pair = number of RT points
    expect_true(res[1, 1, 2L] > 0)
    expect_equal(unname(res[2, 1, 2L]), 0)
})

test_that("compareChromatograms(x, y) with empty objects.", {
    res <- compareChromatograms(c_empty, c_full)
    expect_true(is.array(res))
    expect_equal(dim(res)[1L], 0L)
    res2 <- compareChromatograms(c_full, c_empty)
    expect_true(is.array(res2))
    expect_equal(dim(res2)[2L], 0L)
})

test_that("compareChromatograms(x, y) with custom FUN.", {
    cosine_sim <- function(x, y, ...) {
        sum(x * y, na.rm = TRUE) /
            (sqrt(sum(x^2, na.rm = TRUE)) * sqrt(sum(y^2, na.rm = TRUE)))
    }
    res <- compareChromatograms(c_full[1], c_full[3], FUN = cosine_sim)
    expect_equal(unname(res[1, 1, 1L]), 1)
})

## ---- matchRtime tests ----

test_that("matchRtime returns aligned intensities for overlapping RT", {
    pd_a <- data.frame(rtime = c(1, 2, 3, 4), intensity = c(10, 20, 30, 40))
    pd_b <- data.frame(rtime = c(2, 3, 4, 5), intensity = c(20, 30, 40, 50))
    res <- matchRtime(pd_a, pd_b)
    expect_true(is.list(res))
    expect_equal(length(res$x), length(res$y))
    expect_true(length(res$x) >= 2L)
    ## Overlapping RT is 2-4; both have identical intensities at those points
    expect_equal(res$x, res$y)
})

test_that("matchRtime returns empty vectors for non-overlapping RT", {
    pd_a <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    pd_b <- data.frame(rtime = c(10, 11, 12), intensity = c(40, 50, 60))
    res <- matchRtime(pd_a, pd_b)
    expect_equal(res$x, numeric())
    expect_equal(res$y, numeric())
})

test_that("matchRtime returns empty vectors when fewer than 2 points", {
    pd_short <- data.frame(rtime = 1, intensity = 50)
    pd_ok <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    res <- matchRtime(pd_short, pd_ok)
    expect_equal(res$x, numeric())
    expect_equal(res$y, numeric())
})

test_that("matchRtime interpolates onto union of RT grids", {
    pd_a <- data.frame(rtime = c(1, 3), intensity = c(10, 30))
    pd_b <- data.frame(rtime = c(1, 2, 3), intensity = c(100, 200, 300))
    res <- matchRtime(pd_a, pd_b)
    ## Common grid: 1, 2, 3. pd_a interpolated at rt=2 → 20
    expect_equal(res$x, c(10, 20, 30))
    expect_equal(res$y, c(100, 200, 300))
})

test_that("compareChromatograms with custom MAPFUN.", {
    ## A MAPFUN that returns raw intensities (assumes identical RT)
    identity_map <- function(x, y, ...) {
        list(x = x$intensity, y = y$intensity)
    }
    ## Use chromatograms 1 and 3 which are identical (same number of points)
    res <- compareChromatograms(c_full[c(1, 3)], MAPFUN = identity_map)
    expect_true(is.array(res))
    expect_equal(dim(res), c(2L, 2L, 2L))
    expect_equal(diag(res[, , 1L]), c(1, 1))
    expect_equal(unname(res[1, 2, 1L]), 1)
})
