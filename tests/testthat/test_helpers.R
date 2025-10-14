test_that(".df_combine works as expected", {
    combined_backend <- .df_combine(list(be, be_cd))
    expect_equal(nrow(combined_backend@chromData), nrow(cdata) + nrow(cdata))
    expect_equal(
        length(.peaksData(combined_backend)),
        length(.peaksData(be)) + length(.peaksData(be_cd))
    )

    expect_equal(.df_combine(list(be)), be)
    incompatible_data <- list(
        data.frame(
            rtime = c(10.0, 12.0), intensity = c(200, 150),
            other_col = c("test", "test")
        ),
        data.frame(
            rtime = c(30.1, 31.2), intensity = c(110, 90),
            other_col = c("test", "test")
        )
    )
    be_incompatible <- backendInitialize(be_empty,
        chromData = cdata,
        peaksData = incompatible_data
    )

    expect_error(
        .df_combine(c(be, be_incompatible)),
        "Provided objects have different sets of peak variables."
    )

    setClass("DummyBackend",
        contains = "ChromBackend"
    )
    dm <- new("DummyBackend")

    expect_error(
        .df_combine(c(be, dm)),
        "merge backends of the same type"
    )
})


test_that(".filter_ranges helper function works correctly", {
    query <- data.frame(
        mz = c(112.2, 123.3, 134.4),
        chromIndex = c(1L, 2L, 3L)
    )
    ranges <- numeric()
    expect_error(.filter_ranges(query, ranges, match = "any"), "needs to be")
    ranges <- c(1, 2, 3)
    expect_error(
        .filter_ranges(query, ranges, match = "any"),
        "Length of 'ranges'"
    )

    # All ranges match (full data retained)
    ranges <- c(100, 200, 1, 3)
    res <- .filter_ranges(query, ranges, match = "all")
    expect_equal(res, c(1, 2, 3))
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, c(1, 2, 3))

    # No matches (empty result)
    ranges <- c(500, 600, 4, 5)
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, integer(0))

    # Partially overlapping ranges
    ranges <- c(120, 130, 4, 5)
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, 2)
    res <- .filter_ranges(query, ranges, match = "all")
    expect_equal(res, integer(0))

    # Single row edge case
    query_single <- data.frame(mz = 125.0, chromIndex = 2L)
    ranges <- c(120, 130, 1, 3)
    res <- .filter_ranges(query_single, ranges, match = "all")
    expect_equal(res, 1)
    res <- .filter_ranges(query_single, ranges, match = "any")
    expect_equal(res, 1)

    # Edge case: empty query
    query_empty <- data.frame(mz = numeric(), chromIndex = integer())
    ranges <- c(100, 200, 1, 3)
    res <- .filter_ranges(query_empty, ranges, match = "all")
    expect_equal(res, integer(0))
    res <- .filter_ranges(query_empty, ranges, match = "any")
    expect_equal(res, integer(0))
})

test_that(".check_column_order_and_types works", {
    df_valid <- data.frame(
        rtime = c(1, 2, 3),
        intensity = c(10, 20, 30)
    )
    df_invalid_order <- data.frame(
        intensity = c(10, 20, 30),
        rtime = c(1, 2, 3)
    )
    df_invalid_type <- data.frame(
        rtime = c("a", "b", "c"),
        intensity = c(10, 20, 30)
    )
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.check_column_order_and_types(
        df_valid,
        names(.CORE_PEAKS_VARIABLES),
        .CORE_PEAKS_VARIABLES
    ))
    expect_equal(
        .check_column_order_and_types(
            df_invalid_order,
            names(.CORE_PEAKS_VARIABLES),
            .CORE_PEAKS_VARIABLES
        ),
        "Columns should be in the order 'rtime', 'intensity'."
    )
    expect_equal(
        .check_column_order_and_types(
            df_invalid_type,
            names(.CORE_PEAKS_VARIABLES),
            .CORE_PEAKS_VARIABLES
        ),
        "The peaksData variable(s) rtime have the wrong data type."
    )
    expect_null(.check_column_order_and_types(
        df_empty,
        names(.CORE_PEAKS_VARIABLES),
        .CORE_PEAKS_VARIABLES
    ))
})


test_that(".check_rtime works", {
    df_valid <- data.frame(
        rtime = c(1, 2, 3),
        intensity = c(10, 20, 30)
    )
    df_invalid_na <- data.frame(
        rtime = c(1, NA, 3),
        intensity = c(10, 20, 30)
    )
    df_invalid_increasing <- data.frame(
        rtime = c(3, 2, 1),
        intensity = c(10, 20, 30)
    )
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.check_rtime(df_valid))
    expect_equal(
        .check_rtime(df_invalid_na),
        "'rtime' column contains NA values."
    )
    expect_equal(
        .check_rtime(df_invalid_increasing),
        "'rtime' column is not strictly increasing."
    )
    expect_null(.check_rtime(df_empty))
})

test_that(".validate_entry works", {
    df_valid <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.validate_entry(
        df_valid, 1, names(.CORE_PEAKS_VARIABLES),
        .CORE_PEAKS_VARIABLES
    ))
    expect_null(.validate_entry(
        df_empty, 3,
        names(.CORE_PEAKS_VARIABLES),
        .CORE_PEAKS_VARIABLES
    ))
    expect_equal(
        .validate_entry(
            c(1, 2, 3), 1,
            names(.CORE_PEAKS_VARIABLES),
            .CORE_PEAKS_VARIABLES
        ),
        "Entry 1: all 'peaksData' entries should be of class 'data.frame'"
    )
})

test_that(".run_process_queue ChromBackendMemory work", {
    result <- .run_process_queue(c_empty@backend,
        f = processingChunkFactor(c_empty),
        queue = c_empty@processingQueue
    )
    expect_equal(result, c_empty@backend)

    result <- .run_process_queue(c_full@backend,
        f = processingChunkFactor(c_full),
        queue = c_full@processingQueue
    )
    expect_equal(result, c_full@backend)

    c_queued <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    c_queued <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 200)
    )

    # this test for f = factor() and queue >1
    result <- .run_process_queue(c_queued@backend,
        f = processingChunkFactor(c_queued),
        queue = c_queued@processingQueue
    )
    expect_true(inherits(result, "ChromBackend"))
    peaks_result <- peaksData(result)
    expect_equal(length(peaks_result), length(peaksData(c_full@backend)))
    expect_false(identical(peaks_result, peaksData(c_full@backend)))

    f <- factor(c(1, 1, 2))
    c_queued <- filterPeaksData(c_full,
        variables = c("rtime"),
        ranges = c(12.5, 45.5)
    )
    c_queued <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(100, 200)
    )

    result <- .run_process_queue(c_queued@backend,
        f = f,
        queue = c_queued@processingQueue
    )

    expect_true(inherits(result, "ChromBackend"))
    expect_equal(length(peaksData(result)), length(peaksData(c_full@backend)))

    split_data <- split(c_full@backend, f)
    expect_equal(length(split_data), length(levels(f)))

    f <- factor(c(1, 2))
    tmp <- c_full
    tmp@processingQueue <- list(1)
    expect_error(
        .run_process_queue(tmp@backend,
            f = f,
            queue = tmp@processingQueue
        ),
        "length 'f' has to be equal to the length of 'object'"
    )

    f <- c(1, 2, 3)
    expect_error(
        .run_process_queue(tmp@backend,
            f = f, queue = tmp@processingQueue
        ),
        "f must be a factor"
    )
})

test_that(".run_processing_queue, ChromBackendMzr work", {
    ## without factor and queue == 1
    c_queued <- filterPeaksData(c_mzr,
        variables = c("rtime"),
        ranges = c(12.5, 25.5), keep = FALSE
    )

    result1 <- .run_process_queue(c_queued@backend,
        f = processingChunkFactor(c_queued),
        queue = c_queued@processingQueue
    )

    expect_true(result1@inMemory)
    expect_false(c_queued@backend@inMemory)

    expect_true(inherits(result1, "ChromBackendMzR"))
    expect_false(identical(
        lengths(rtime(result1)),
        lengths(rtime(c_mzr@backend))
    ))
    expect_equal(length(peaksData(result1)), length(peaksData(c_mzr@backend)))

    ## with levels(factor) > 1 and queue == 1
    processingChunkSize(c_queued) <- 100
    f <- processingChunkFactor(c_queued) # > 1
    result2 <- .run_process_queue(c_queued@backend,
        f = f,
        queue = c_queued@processingQueue
    )
    expect_true(inherits(result2, "ChromBackendMzR"))
    expect_equal(length(peaksData(result2)), length(peaksData(c_mzr@backend)))
    expect_false(identical(
        lengths(rtime(result2)),
        lengths(rtime(c_mzr@backend))
    ))
    expect_true(result2@inMemory)
    expect_false(c_queued@backend@inMemory)
    expect_identical(result1, result2)


    ## without factor and queue > 1
    c_queued <- filterPeaksData(c_queued,
        variables = c("intensity"),
        ranges = c(45, 50)
    )
    result3 <- .run_process_queue(c_queued@backend,
        f = factor(),
        queue = c_queued@processingQueue
    )

    expect_true(inherits(result3, "ChromBackendMzR"))
    expect_equal(length(peaksData(result3)), length(peaksData(c_mzr@backend)))
    expect_false(identical(
        lengths(rtime(result3)),
        lengths(rtime(c_mzr@backend))
    ))
    expect_true(result3@inMemory)
    expect_false(c_queued@backend@inMemory)


    ## with factor and queue > 1
    f <- processingChunkFactor(c_queued)
    result4 <- .run_process_queue(c_queued@backend,
        f = f,
        queue = c_queued@processingQueue
    )
    expect_true(inherits(result4, "ChromBackendMzR"))
    expect_equal(
        length(peaksData(result4)),
        length(peaksData(c_mzr@backend))
    )
    expect_false(identical(
        lengths(rtime(result4)),
        lengths(rtime(c_mzr@backend))
    ))
    expect_true(result4@inMemory)
    expect_false(c_queued@backend@inMemory)
    expect_identical(result3, result4)
})

test_that(".valid_processing_queue works correctly", {
    expect_null(.valid_processing_queue(list()))
    valid_queue <- list(new("ProcessingStep"))
    expect_null(.valid_processing_queue(valid_queue))

    invalid_queue <- list("not_a_processing_step")
    expect_error(
        .valid_processing_queue(invalid_queue),
        "'processingQueue' should only contain ProcessingStep objects."
    )
})

test_that("ensure_rt_mz_columns correctly handles mz and rt columns", {
    spectra <- s
    spectra_f <- factor(
        do.call(
            paste,
            c(as.list(Spectra::spectraData(s)[, c("msLevel", "dataOrigin")]),
              sep = "_")))
    chrom_data <- data.frame(msLevel = c(1,2,3))
    chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
    expect_equal(chrom_data$mzMin, c(-Inf, -Inf, -Inf))
    expect_equal(chrom_data$mzMax, c(Inf, Inf, Inf))

    chrom_data <- data.frame(mzMin = c(100))
    expect_error(.ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
                 "must be present if one is provided.")

    chrom_data <- data.frame(mzMax = c(200))
    expect_error(.ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
                 "must be present if one is provided.")
    chrom_data <- data.frame(msLevel = c(1,2,3))
    chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
    s_plit <- split(spectra, spectra_f)
    expect_equal(chrom_data$rtMin[[1]], min(s_plit[[1]]$rtime, na.rm = TRUE))
    expect_equal(chrom_data$rtMax[[1]], max(s_plit[[1]]$rtime, na.rm = TRUE))

    chrom_data <- data.frame(rtMin = c(10))
    expect_error(.ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
                 " must be present if one is provided.")
    chrom_data <- data.frame(rtMax = c(50))
    expect_error(.ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
                 "must be present if one is provided.")

    chrom_data <- data.frame(mzMin = c(100), mzMax = c(200),
                             rtMin = c(10), rtMax = c(50))
    chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
    expect_equal(chrom_data$mzMin, 100)
    expect_equal(chrom_data$mzMax, 200)
    expect_equal(chrom_data$rtMin, 10)
    expect_equal(chrom_data$rtMax, 50)
})

test_that(".validate_chromExtract_input works correctly", {
    cdata <- data.frame(
        msLevel = c(1L, 1L, 1L),
        dataOrigin = c("A", "B", "C")
    )
    be <- backendInitialize(new("ChromBackendMemory"), chromData = cdata)

    peak_tbl <- data.frame(
        rtMin = c(1, 2, 3),
        rtMax = c(5, 6, 7),
        msLevel = c(1L, 1L, 1L),
        dataOrigin = c("A", "B", "C")
    )

    # should pass
    expect_silent(
        .validate_chromExtract_input(be, peak_tbl,
                                     by = c("msLevel", "dataOrigin"))
    )

    # missing required column
    bad_tbl <- peak_tbl[, !names(peak_tbl) %in% "rtMax", drop = FALSE]
    expect_error(
        .validate_chromExtract_input(be, bad_tbl,
                                     by = c("msLevel", "dataOrigin")),
        "must contain columns"
    )

    # NA in rtMin
    bad_tbl2 <- peak_tbl
    bad_tbl2$rtMin[1] <- NA
    expect_error(
        .validate_chromExtract_input(be, bad_tbl2,
                                     by = c("msLevel", "dataOrigin")),
        "cannot contain NA"
    )

    # missing 'by' columns in chromData
    bad_be <- backendInitialize(new("ChromBackendMemory"),
                                chromData = cdata[, "msLevel", drop = FALSE])
    expect_error(
        .validate_chromExtract_input(bad_be, peak_tbl,
                                     by = c("msLevel", "dataOrigin")),
        "must be present"
    )

    ## unique
    expect_error(
        .validate_chromExtract_input(bad_be, peak_tbl,
                                     by = "msLevel"),
        "must uniquely identify rows"
    )

})

test_that(".match_chromdata_peaktable aligns correctly", {
    tmp_cdata <- data.frame(
        msLevel = c(1L, 1L, 2L),
        dataOrigin = c("A", "B", "A")
    )
    tmp <- backendInitialize(new("ChromBackendMemory"), chromData = tmp_cdata)

    peak_tbl <- data.frame(
        msLevel = c(1L, 2L),
        dataOrigin = c("A", "A"),
        rtMin = c(1, 2),
        rtMax = c(5, 6)
    )

    matched <- .match_chromdata_peaktable(tmp, peak_tbl,
                                          by = c("msLevel", "dataOrigin"))

    # Expect a subset of object
    expect_s4_class(matched$object, "ChromBackendMemory")
    expect_equal(length(matched$chrom_keys), nrow(.chromData(matched$object)))

    # Check factor levels alignment
    expect_true(all(levels(matched$peak_keys) %in% levels(matched$chrom_keys)))

    # missing key should error
    bad_tbl <- data.frame(
        msLevel = 3L, dataOrigin = "Z", rtMin = 1, rtMax = 2
    )
    expect_error(
        .match_chromdata_peaktable(tmp, bad_tbl,
                                   by = c("msLevel", "dataOrigin")),
        "do not exist"
    )
})

test_that(".check_overl_columns warns correctly", {
    tmp_cdata <- data.frame(
        msLevel = 1L,
        dataOrigin = "X",
        mz = 100, extracol = "info"
    )
    tmp <- backendInitialize(new("ChromBackendMemory"), chromData = tmp_cdata)

    peak_tbl <- data.frame(
        rtMin = 1, rtMax = 2, mzMin = 99, mzMax = 101,
        msLevel = 1L, dataOrigin = "X", mz = 123, extracol = "test"
    )

    req_cols <- c("rtMin", "rtMax", "mzMin", "mzMax", "msLevel", "dataOrigin")

    expect_warning(
        overl <- .check_overl_columns(tmp, peak_tbl, req_cols),
        "already exist"
    )

    # overlapping should include "mz"
    expect_true(all(c("mz", "extracol") %in% names(peak_tbl)[overl]))

})

test_that(".impute() works correctly and without warnings", {
    # Base signal with gaps
    x <- c(1:5, NA, 7:10, NA, 12:15, rep(NA, 2), 18:20)

    ## linear
    expect_silent({
        res_lin <- .impute(x, method = "linear")
    })
    expect_false(anyNA(res_lin))
    expect_true(all(diff(res_lin) > 0))  # still increasing

    ## spline
    expect_silent({
        res_spl <- .impute(x, method = "spline")
    })
    expect_false(anyNA(res_spl))
    expect_equal(length(res_spl), length(x))

    ## Gaussian
    expect_silent(
        res_gauss <- .impute(x, method = "gaussian", window = 2, sd = 1)
    )
    expect_false(anyNA(res_gauss))
    expect_equal(length(res_gauss), length(x))

    ## loess
    expect_warning({
        res_loess <- .impute(x, method = "loess", span = 0.3)
        "could not fill all NAs"
    })
    expect_false(anyNA(res_loess))
    expect_equal(length(res_loess), length(x))

    ## Consecutive NAs
    x_na <- c(1, 2, NA, NA, 5, 6, 7, 8, NA, NA, 11, 12)
    expect_silent({
        res_consec <- .impute(x_na, method = "linear")
    })
    expect_false(anyNA(res_consec))
    expect_equal(length(res_consec), length(x_na))

    ## No NA
    x_nomiss <- 1:10
    expect_silent({
        res_nomiss <- .impute(x_nomiss, method = "spline")
    })
    expect_identical(res_nomiss, x_nomiss)

    ## all Nas returns NA
    x_allna <- rep(NA_real_, 8)
    expect_silent({
        res_allna <- .impute(x_allna, method = "gaussian")
    })
    expect_true(all(is.na(res_allna)))
})



