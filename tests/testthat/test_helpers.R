test_that(".df_combine works as expected", {
    combined_backend <- .df_combine(list(be, be_cd))
    expect_equal(nrow(combined_backend@chromData), nrow(cdata) + nrow(cdata))
    expect_equal(length(combined_backend@peaksData),
                 length(be@peaksData) + length(be_cd@peaksData))

    expect_equal(.df_combine(list(be)), be)
    incompatible_data <- list(
        data.frame(rtime = c(10.0, 12.0), intensity = c(200, 150),
                   other_col= c("test", "test")),
        data.frame(rtime = c(30.1, 31.2), intensity = c(110, 90),
                   other_col= c("test", "test"))
    )
    be_incompatible <- backendInitialize(be_empty, chromData = cdata,
                                         peaksData = incompatible_data)

    expect_error(.df_combine(c(be, be_incompatible)),
                 "Provided objects have different sets of peak variables.")

    setClass("DummyBackend",
             contains = "ChromBackend")
    dm <- new("DummyBackend")

    expect_error(.df_combine(c(be, dm)),
                 "merge backends of the same type")
})


test_that(".filter_ranges helper function works correctly", {
    query <- data.frame(mz = c(112.2, 123.3, 134.4), chromIndex = c(1L, 2L, 3L))
    ranges <- numeric()
    expect_error(.filter_ranges(query, ranges, match = "any"), "needs to be")
    ranges <- c(1, 2, 3)
    expect_error(.filter_ranges(query, ranges, match = "any"), "Length of 'ranges'")

    # All ranges match (full data retained)
    ranges <- c(100, 200, 1, 3)
    res <- .filter_ranges(query, ranges, match = "all")
    expect_equal(res, c(1,2,3))
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, c(1,2,3))

    # No matches (empty result)
    ranges <- c(500, 600, 4,5)
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, integer(0))

    # Partially overlapping ranges
    ranges <- c(120, 130, 4, 5)
    res <- .filter_ranges(query, ranges, match = "any")
    expect_equal(res, 2)
    res <- .filter_ranges(query, ranges, match = "all")
    expect_equal(res, integer(0))
})


test_that(".valid_chrom_backend_data_storage works", {
    expect_match(.valid_chrom_backend_data_storage(c("a", NA)), "not allowed")
    expect_null(.valid_chrom_backend_data_storage(character()))
    expect_null(.valid_chrom_backend_data_storage("b"))
})

test_that(".check_column_order_and_types works", {
    df_valid <- data.frame(rtime = c(1, 2, 3),
                           intensity = c(10, 20, 30))
    df_invalid_order <- data.frame(intensity = c(10, 20, 30),
                                   rtime = c(1, 2, 3))
    df_invalid_type <- data.frame(rtime = c("a", "b", "c"),
                                  intensity = c(10, 20, 30))
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.check_column_order_and_types(df_valid,
                                              names(.CORE_PEAKS_VARIABLES),
                                              .CORE_PEAKS_VARIABLES))
    expect_equal(.check_column_order_and_types(df_invalid_order,
                                               names(.CORE_PEAKS_VARIABLES),
                                               .CORE_PEAKS_VARIABLES),
                 "Columns should be in the order 'rtime', 'intensity'.")
    expect_equal(.check_column_order_and_types(df_invalid_type,
                                               names(.CORE_PEAKS_VARIABLES),
                                               .CORE_PEAKS_VARIABLES),
                 "The peaksData variable(s) rtime have the wrong data type.")
    expect_null(.check_column_order_and_types(df_empty,
                                              names(.CORE_PEAKS_VARIABLES),
                                              .CORE_PEAKS_VARIABLES))
})


test_that(".check_rtime works", {
    df_valid <- data.frame(rtime = c(1, 2, 3),
                           intensity = c(10, 20, 30))
    df_invalid_na <- data.frame(rtime = c(1, NA, 3),
                                intensity = c(10, 20, 30))
    df_invalid_increasing <- data.frame(rtime = c(3, 2, 1),
                                        intensity = c(10, 20, 30))
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.check_rtime(df_valid))
    expect_equal(.check_rtime(df_invalid_na),
                 "'rtime' column contains NA values.")
    expect_equal(.check_rtime(df_invalid_increasing),
                 "'rtime' column is not strictly increasing.")
    expect_null(.check_rtime(df_empty))
})

test_that(".validate_entry works", {
    df_valid <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    df_empty <- data.frame(rtime = numeric(), intensity = numeric())

    expect_null(.validate_entry(df_valid, 1, names(.CORE_PEAKS_VARIABLES),
                                .CORE_PEAKS_VARIABLES))
    expect_null(.validate_entry(df_empty, 3,
                                names(.CORE_PEAKS_VARIABLES),
                                .CORE_PEAKS_VARIABLES))
    expect_equal(.validate_entry(c(1,2,3), 1,
                                 names(.CORE_PEAKS_VARIABLES),
                                 .CORE_PEAKS_VARIABLES),
                 "Entry 1: all 'peaksData' entries should be of class 'data.frame'")
})

test_that(".run_process_queue handles empty queue correctly", {
    result <- .run_process_queue(c_empty@backend,
                                 f = processingChunkFactor(c_empty),
                                 queue = c_empty@processingQueue)
    expect_equal(result, c_empty@backend)

    result <- .run_process_queue(c_full@backend,
                                 f = processingChunkFactor(c_full),
                                 queue = c_full@processingQueue)
    expect_equal(result, c_full@backend)

    c_queued <- filterPeaksData(c_full, variables = c("rtime"), ranges = c(12.5, 45.5))
    c_queued <- filterPeaksData(c_queued, variables = c("intensity"), ranges = c(100, 200))

    result <- .run_process_queue(c_queued@backend,
                                 f = processingChunkFactor(c_queued), # this test for f = factor() and queue >1
                                 queue = c_queued@processingQueue)
    expect_true(inherits(result, "ChromBackend"))
    peaks_result <- peaksData(result)
    expect_equal(length(peaks_result), length(peaksData(c_full@backend)))
    expect_false(identical(peaks_result, peaksData(c_full@backend)))

    f <- factor(c(1, 1, 2))
    c_queued <- filterPeaksData(c_full, variables = c("rtime"), ranges = c(12.5, 45.5))
    c_queued <- filterPeaksData(c_queued, variables = c("intensity"), ranges = c(100, 200))

    result <- .run_process_queue(c_queued@backend, f = f,
                                 queue = c_queued@processingQueue)

    expect_true(inherits(result, "ChromBackend"))
    expect_equal(length(peaksData(result)), length(peaksData(c_full@backend)))

    split_data <- split(c_full@backend, f)
    expect_equal(length(split_data), length(levels(f)))

    f <- factor(c(1, 2))
    tmp <- c_full
    tmp@processingQueue <- list(1)
    expect_error(.run_process_queue(tmp@backend, f = f, queue <- tmp@processingQueue),
                 "length 'f' has to be equal to the length of 'object'")

    f <- c(1, 2, 3)
    expect_error(.run_process_queue(tmp@backend,
                                    f = f, queue = tmp@processingQueue),
                 "f must be a factor")

    f <- factor(c(1, 1, 2))
    c_empty <- Chromatograms()
    result <- .run_process_queue(c_empty@backend, f = f,
                                 queue = c_empty@processingQueue)
    expect_equal(result, c_empty@backend)
})


