test_that(".valid_chrom_data_required_columns works", {
    df <- DataFrame()
    expect_null(.valid_chrom_data_required_columns(df))
    df <- DataFrame(msLevel = 1L)
    expect_match(.valid_chrom_data_required_columns(df),
                 "Required column")
    df$dataStorage <- "some"
    expect_null(.valid_chrom_data_required_columns(df))
})

test_that(".valid_column_datatype works", {
    expect_null(.valid_column_datatype(DataFrame(a = 4)))
    expect_match(.valid_column_datatype(
        DataFrame(msLevel = "a", mz = TRUE, intensity = list(1:4))),
        "type: msLevel, mz.")
    expect_null(.valid_column_datatype(
        DataFrame(msLevel = 1L, mz = 2.4)))
})

test_that(".valid_rtime_column works", {
    df <- DataFrame(msLevel = c(1L, 1L))
    df$rtime <- IRanges::NumericList(1:3, 1:12)
    expect_null(.valid_rtime_column(df))
    expect_null(.valid_rtime_column(DataFrame(msLevel = 4L)))
    df$rtime <- list(1:3, letters[1:4])
    expect_match(.valid_rtime_column(df),
                 "rtime column should be of type NumericList")
    df$rtime <- IRanges::NumericList(1:4, c(3, 2, 5))
    expect_match(.valid_rtime_column(df),
                 "sorted increasingly")
})

test_that(".valid_intensity_column works", {
    df <- DataFrame(msLevel = c(1L, 1L))
    df$intensity <- IRanges::NumericList(4, 2:5)
    expect_null(.valid_intensity_column(df))
    df$intensity <- list("g", TRUE)
    expect_match(.valid_intensity_column(df),
                 "intensity column should be of type NumericList")
})

test_that(".valid_intensity_rtime_columns works", {
    be <- ChromBackendDataFrame()
    expect_null(.valid_intensity_rtime_columns(be@chromData))
    be <- backendInitialize(be, DataFrame(fromFile = c(1L, 1L)))
    expect_null(.valid_intensity_rtime_columns(be@chromData))
    be@chromData$rtime <- list(1:3, 1:2)
    be@chromData$intensity <- list(1:3, 2)
    expect_match(.valid_intensity_rtime_columns(be@chromData),
                 "Length of rtime and intensity")
    be@chromData$intensity <- list(1:3, 1:2)
    expect_null(.valid_intensity_rtime_columns(be@chromData))
})

test_that(".get_column works", {
    df <- DataFrame(msLevel = c(NA_integer_, NA_integer_),
                    fromFile = rep(1L, 2), other_col = rep("a", 2))
    res <- .get_column(df, column = "fromFile")
    expect_equal(res, c(1L, 1L))
    res <- .get_column(df, column = "other_col")
    expect_equal(res, c("a", "a"))
    res <- .get_column(df, column = "msLevel")
    expect_equal(res, c(NA_integer_, NA_integer_))
    res <- .get_column(df, column = "productMz")
    expect_equal(res, c(NA_real_, NA_real_))
    expect_error(.get_column(df, column = "a"), "not available")
    df <- DataFrame()
    res <- .get_column(df, column = "msLevel")
    expect_equal(res, integer())
})

test_that(".combine_chrom_backend_data_frame works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), fromFile = 1L,
                    mz = as.numeric(1:3))
    df2 <- DataFrame(msLevel = c(2L, 1L), fromFile = 1L,
                     mz = c(4.1, 5.2), scanIndex = 1:2)
    df3 <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L,
                     precScanNum = 1L, other_col = "z")
    be <- backendInitialize(ChromBackendDataFrame(), df)
    be2 <- backendInitialize(ChromBackendDataFrame(), df2)
    be3 <- backendInitialize(ChromBackendDataFrame(), df3)

    expect_equal(.combine_chrom_backend_data_frame(list(be)), be)
    expect_error(backendMerge(list(be, 4)), "backends of the same type")

    res <- .combine_chrom_backend_data_frame(list(be, be2, be3))
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_identical(res@chromData$dataStorage, rep("<memory>", 7))
    expect_identical(dataStorage(res), rep("<memory>", 7))
    expect_identical(msLevel(res), c(1L, 2L, 2L, 2L, 1L, 1L, 2L))
    expect_identical(mz(res), c(1:3, 4.1, 5.2, NA, NA))
    expect_identical(res@chromData$other_col,
                     c(rep(NA_character_, 5), "z", "z"))
    expect_true(is(be3@chromData$precScanNum, "integer"))
    expect_true(is(res@chromData$precScanNum, "integer"))

    ## One backend with and one without rtime
    df2$rtime <- list(c(1.1, 1.2), c(1.1, 1.2))
    df2$intensity <- list(c(12.4, 3), c(123.4, 1))
    be2 <- backendInitialize(ChromBackendDataFrame(), df2)
    res <- .combine_chrom_backend_data_frame(list(be, be2, be3))
    expect_equal(lengths(rtime(res)), c(0, 0, 0, 2, 2, 0, 0))

    ## With different dataStorage
    be$dataStorage <- c("a", "a", "a")
    be3$dataStorage <- c("z", "b")

    res <- .combine_chrom_backend_data_frame(list(be, be2, be3))
    expect_identical(res$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(res@chromData$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(mz(res), c(1:3, 4.1, 5.2, NA, NA))
})
