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


test_that(".check_rtime wprks", {
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

