test_df <- DataFrame(msLevel = c(1L, 2L, 2L), mz = as.numeric(4:6))
test_df$rtime <- list(c(1.1, 1.3, 1.5), c(4.1, 5.1), c(1.6, 1.7, 1.8, 1.9))
test_df$intensity <- list(c(45.1, 34, 12), c(234.4, 1333), c(42.1, 34.2, 65, 6))

test_that("show,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    show(be)
    df <- DataFrame(rt = c(1.2, 1.3))
    be <- backendInitialize(be, df)
    show(be)
})

test_that("backendInitialize,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_true(validObject(be))
    be <- backendInitialize(be)
    expect_true(validObject(be))
    be <- backendInitialize(be, chromData = DataFrame(msLevel = 2L))
    expect_true(validObject(be))
    expect_equal(dataStorage(be), "<memory>")

    df <- test_df
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@chromData$rtime, "NumericList"))
    expect_true(is(be@chromData$intensity, "NumericList"))

    df$rtime <- SimpleList(df$rtime)
    df$intensity <- SimpleList(df$intensity)
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@chromData$rtime, "NumericList"))
    expect_true(is(be@chromData$intensity, "NumericList"))
    expect_identical(be@chromData$dataStorage, rep("<memory>", 3))
    expect_identical(be$dataStorage, rep("<memory>", 3))

    df$rtime <- NumericList(df$rtime)
    df$intensity <- NumericList(df$intensity)
    be <- backendInitialize(be, df)
    expect_true(validObject(be))
    expect_true(is(be@chromData$rtime, "NumericList"))
    expect_true(is(be@chromData$intensity, "NumericList"))
})

test_that("backendMerge,ChromBackendDataFrame works", {
    df <- DataFrame(msLevel = c(1L, 2L, 2L), fromFile = 1L,
                    mz = as.numeric(1:3))
    df2 <- DataFrame(msLevel = c(2L, 1L), fromFile = 1L,
                     mz = c(4.1, 5.2), scanIndex = 1:2)
    df3 <- DataFrame(msLevel = c(1L, 2L), fromFile = 1L,
                     precScanNum = 1L, other_col = "z")
    be <- backendInitialize(ChromBackendDataFrame(), df)
    be2 <- backendInitialize(ChromBackendDataFrame(), df2)
    be3 <- backendInitialize(ChromBackendDataFrame(), df3)

    expect_equal(backendMerge(be), be)
    expect_error(backendMerge(be, 4), "backends of the same type")

    res <- backendMerge(be, be2, be3)
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
    res <- backendMerge(be, be2, be3)
    expect_identical(lengths(rtime(res)), c(0L, 0L, 0L, 2L, 2L, 0L, 0L))
    expect_identical(lengths(rtime(res)), lengths(intensity(res)))

    ## With different dataStorage
    be$dataStorage <- c("a", "a", "a")
    be3$dataStorage <- c("z", "b")

    res <- backendMerge(be, be2, be3)
    expect_identical(res$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(res@chromData$dataStorage,
                     c("a", "a", "a", "<memory>", "<memory>", "z", "b"))
    expect_identical(mz(res), c(1:3, 4.1, 5.2, NA, NA))
})

test_that("chromData, chromData<-, ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    res <- chromData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_equal(colnames(res), names(.CHROMATOGRAMS_DATA_COLUMNS))

    df <- DataFrame(mz = c(1.2, 1.4), a = "a", b = "b")
    be <- backendInitialize(be, chromData = df)

    res <- chromData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(all(names(.CHROMATOGRAMS_DATA_COLUMNS) %in% colnames(res)))
    expect_equal(res$a, c("a", "a"))
    expect_equal(res$b, c("b", "b"))

    res <- chromData(be, "msLevel")
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "msLevel")
    expect_equal(res$msLevel, c(NA_integer_, NA_integer_))

    res <- chromData(be, c("rtime"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), "rtime")
    expect_equal(res$rtime, NumericList(numeric(), numeric(), compress = FALSE))

    res <- chromData(be, c("a", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), c("a", "intensity"))
    expect_equal(res$intensity, NumericList(numeric(), numeric(),
                                            compress = FALSE))
    expect_equal(res$a, c("a", "a"))

    chromData(be) <- DataFrame(mzLevel = c(3L, 4L),
                               mz = c(1.2, 1.4), other_col = "b")
    expect_identical(be$mz, c(1.2, 1.4))
    expect_true(any(chromVariables(be) == "other_col"))
    expect_identical(chromData(be, "other_col")[, 1], c("b", "b"))

    expect_error(chromData(be) <- DataFrame(msLevel = 1:3),
                 "with 2 rows")
})

test_that("chromIndex works", {
    df <- test_df
    df$chromIndex <- c("a", "b", "c")
    expect_error(backendInitialize(ChromBackendDataFrame(), df),
                 "wrong data")
    df$chromIndex <- 1:3
    be <- backendInitialize(ChromBackendDataFrame(), df)
    expect_identical(chromIndex(be), 1:3)
})

test_that("chromNames, chromNames<-,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_null(chromNames(be))
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, chromData = df)
    expect_null(chromNames(be))
    df <- DataFrame(msLevel = c(1L, 2L))
    rownames(df) <- c("sp_1", "sp_2")
    be <- backendInitialize(be, chromData = df)
    expect_equal(chromNames(be), c("sp_1", "sp_2"))
    expect_error(chromNames(be) <- "a", "rownames length")
    chromNames(be) <- c("a", "b")
    expect_equal(chromNames(be), c("a", "b"))
})

test_that("chromVariables,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(chromVariables(be), names(.CHROMATOGRAMS_DATA_COLUMNS))
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, chromData = df)
    expect_equal(chromVariables(be), names(.CHROMATOGRAMS_DATA_COLUMNS))
    df$other_column <- 3
    be <- backendInitialize(be, chromData = df)
    expect_equal(chromVariables(be), c(names(.CHROMATOGRAMS_DATA_COLUMNS),
                                         "other_column"))
})

test_that("collisionEnergy,collisionEnergy<-,ChromBackendDataFrame work", {
    be <- ChromBackendDataFrame()
    expect_true(is.numeric(collisionEnergy(be)))
    be <- backendInitialize(be, test_df)
    expect_identical(collisionEnergy(be), rep(NA_real_, length(be)))

    df <- test_df
    df$collisionEnergy <- c(12.2, NA_real_, 34.4)
    be <- backendInitialize(be, df)
    expect_identical(collisionEnergy(be), c(12.2, NA_real_, 34.4))

    collisionEnergy(be) <- c(34.1, 23.5, 1.2)
    expect_identical(collisionEnergy(be), c(34.1, 23.5, 1.2))
})

test_that("dataOrigin,dataOrigin<-,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(dataOrigin(be), character())
    be <- backendInitialize(be, chromData = DataFrame(msLevel = c(1L, 2L)))
    expect_identical(dataOrigin(be), rep(NA_character_, 2))
    expect_error(dataOrigin(be) <- "a", "of length 2")
    dataOrigin(be) <- c("b", "a")
    expect_identical(dataOrigin(be), c("b", "a"))
})

test_that("dataStorage,dataStorage<-,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(dataStorage(be), character())
    be <- backendInitialize(be, chromData = DataFrame(msLevel = c(1L, 2L)))
    expect_identical(dataStorage(be), rep("<memory>", 2))
    dataStorage(be) <- c("a", "b")
    expect_identical(dataStorage(be), c("a", "b"))
    expect_error(dataStorage(be) <- c("a", NA), "not allowed")
})

test_that("intensity,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(intensity(be), NumericList(compress = FALSE))
    be <- backendInitialize(be, chromData = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(intensity(be), NumericList(numeric(), numeric(),
                                            compress = FALSE))
    df <- DataFrame(msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$rtime <- list(1:4, 1:2)
    be <- backendInitialize(be, chromData = df)
    expect_equal(intensity(be), NumericList(1:4, c(2.1, 3.4), compress = FALSE))
})

test_that("intensity<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(), chromData = test_df)

    new_ints <- lapply(test_df$intensity, function(z) z / 2)
    intensity(be) <- new_ints
    expect_identical(intensity(be), NumericList(new_ints, compress = FALSE))

    expect_error(intensity(be) <- 3, "has to be a list")
    expect_error(intensity(be) <- list(3, 2), "match the length")
    expect_error(intensity(be) <- list(3, 2, 4), "number of data pairs")
})

test_that("isEmpty,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(isEmpty(be), logical())
    df <- DataFrame(msLevel = c(1L, 2L))
    be <- backendInitialize(be, chromData = df)
    expect_equal(isEmpty(be), c(TRUE, TRUE))
    df$intensity <- list(1:2, 1:5)
    df$rtime <- list(1:2, 1:5)
    be <- backendInitialize(be, chromData = df)
    expect_equal(isEmpty(be), c(FALSE, FALSE))
})

test_that("length,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(length(be), 0)
    be <- new("ChromBackendDataFrame", chromData = DataFrame(a = 1:3,
                                                             dataStorage = "a"))
    expect_equal(length(be), 3)
})

test_that("msLevel,msLevel<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(msLevel(be), c(1, 2, 1))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(scanIndex = 1:4))
    expect_equal(msLevel(be), rep(NA_integer_, 4))
    msLevel(be) <- c(2, 4, 2, 1)
    expect_identical(msLevel(be), c(2L, 4L, 2L, 1L))

    expect_error(msLevel(be) <- c(2, 3, 1), "of length 4")
    expect_error(msLevel(be) <- c("a", "b", "d", "d"), "of length 4")
})

test_that("mz,mz<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(mz = c(1.2, 1.4)))
    expect_equal(mz(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(mz(be), rep(NA_real_, 3))
    mz(be) <- c(1.2, 1.4, 1.5)
    expect_identical(mz(be), c(1.2, 1.4, 1.5))

    expect_error(mz(be) <- c(2, 3), "of length 3")
    expect_error(mz(be) <- c("a", "b", "d"), "of length 3")
})

test_that("mzMax,mzMax<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(mzMax = c(1.2, 1.4)))
    expect_equal(mzMax(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(mzMax(be), rep(NA_real_, 3))
    mz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(mzMax(be), c(3.1, 4.3, 1.2))
    mzMax(be) <- c(1.2, 1.4, 1.5)
    expect_identical(mzMax(be), c(1.2, 1.4, 1.5))

    expect_error(mzMax(be) <- c(2, 3), "of length 3")
    expect_error(mzMax(be) <- c("a", "b", "d"), "of length 3")
})

test_that("mzMin,mzMin<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(mzMin = c(1.2, 1.4)))
    expect_equal(mzMin(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(mzMin(be), rep(NA_real_, 3))
    mz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(mzMin(be), c(3.1, 4.3, 1.2))
    mzMin(be) <- c(1.2, 1.4, 1.5)
    expect_identical(mzMin(be), c(1.2, 1.4, 1.5))

    expect_error(mzMin(be) <- c(2, 3), "of length 3")
    expect_error(mzMin(be) <- c("a", "b", "d"), "of length 3")
})

test_that("as.list,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(as.list(be), list())
    df <- DataFrame(msLevel = c(1L, 1L))
    be <- backendInitialize(be, chromData = df)
    expect_equal(as.list(be),
                 list(cbind(rtime = numeric(), intensity = numeric()),
                      cbind(rtime = numeric(), intensity = numeric())))
    df$rtime <- list(1:3, c(2.1))
    df$intensity <- list(1:3, 4)
    be <- backendInitialize(be, chromData = df)
    expect_equal(as.list(be), list(cbind(rtime = 1:3, intensity = 1:3),
                                   cbind(rtime = 2.1, intensity = 4)))
})

test_that("lenghts,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(), test_df)
    expect_equal(lengths(be), c(3, 2, 4))

    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(lengths(be), c(0, 0, 0))
})

test_that("precursorMz,precursorMz<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(precursorMz = c(1.2, 1.4)))
    expect_equal(precursorMz(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(precursorMz(be), rep(NA_real_, 3))
    precursorMz(be) <- c(1.2, 1.4, 1.5)
    expect_identical(precursorMz(be), c(1.2, 1.4, 1.5))

    expect_error(precursorMz(be) <- c(2, 3), "of length 3")
    expect_error(precursorMz(be) <- c("a", "b", "d"), "of length 3")
})

test_that("precursorMzMax,precursorMzMax<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(precursorMzMax = c(1.2, 1.4)))
    expect_equal(precursorMzMax(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(precursorMzMax(be), rep(NA_real_, 3))
    precursorMz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(precursorMzMax(be), c(3.1, 4.3, 1.2))
    precursorMzMax(be) <- c(1.2, 1.4, 1.5)
    expect_identical(precursorMzMax(be), c(1.2, 1.4, 1.5))

    expect_error(precursorMzMax(be) <- c(2, 3), "of length 3")
    expect_error(precursorMzMax(be) <- c("a", "b", "d"), "of length 3")
})

test_that("precursorMzMin,precursorMzMin<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(precursorMzMin = c(1.2, 1.4)))
    expect_equal(precursorMzMin(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(precursorMzMin(be), rep(NA_real_, 3))
    precursorMz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(precursorMzMin(be), c(3.1, 4.3, 1.2))
    precursorMzMin(be) <- c(1.2, 1.4, 1.5)
    expect_identical(precursorMzMin(be), c(1.2, 1.4, 1.5))

    expect_error(precursorMzMin(be) <- c(2, 3), "of length 3")
    expect_error(precursorMzMin(be) <- c("a", "b", "d"), "of length 3")
})

test_that("productMz,productMz<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(productMz = c(1.2, 1.4)))
    expect_equal(productMz(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(productMz(be), rep(NA_real_, 3))
    productMz(be) <- c(1.2, 1.4, 1.5)
    expect_identical(productMz(be), c(1.2, 1.4, 1.5))

    expect_error(productMz(be) <- c(2, 3), "of length 3")
    expect_error(productMz(be) <- c("a", "b", "d"), "of length 3")
})

test_that("productMzMax,productMzMax<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(productMzMax = c(1.2, 1.4)))
    expect_equal(productMzMax(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(productMzMax(be), rep(NA_real_, 3))
    productMz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(productMzMax(be), c(3.1, 4.3, 1.2))
    productMzMax(be) <- c(1.2, 1.4, 1.5)
    expect_identical(productMzMax(be), c(1.2, 1.4, 1.5))

    expect_error(productMzMax(be) <- c(2, 3), "of length 3")
    expect_error(productMzMax(be) <- c("a", "b", "d"), "of length 3")
})

test_that("productMzMin,productMzMin<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(productMzMin = c(1.2, 1.4)))
    expect_equal(productMzMin(be), c(1.2, 1.4))
    be <- backendInitialize(ChromBackendDataFrame(),
                            DataFrame(msLevel = c(1L, 2L, 1L)))
    expect_equal(productMzMin(be), rep(NA_real_, 3))
    productMz(be) <- c(3.1, 4.3, 1.2)
    expect_equal(productMzMin(be), c(3.1, 4.3, 1.2))
    productMzMin(be) <- c(1.2, 1.4, 1.5)
    expect_identical(productMzMin(be), c(1.2, 1.4, 1.5))

    expect_error(productMzMin(be) <- c(2, 3), "of length 3")
    expect_error(productMzMin(be) <- c("a", "b", "d"), "of length 3")
})

test_that("rtime,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(rtime(be), NumericList(compress = FALSE))
    be <- backendInitialize(be, chromData = DataFrame(msLevel = c(1L, 2L)))
    expect_equal(rtime(be), NumericList(numeric(), numeric(),
                                        compress = FALSE))
    df <- DataFrame(msLevel = c(1L, 2L))
    df$intensity <- list(1:4, c(2.1, 3.4))
    df$rtime <- list(1:4, 1:2)
    be <- backendInitialize(be, chromData = df)
    expect_equal(rtime(be), NumericList(1:4, 1:2, compress = FALSE))
})

test_that("rtime<-,ChromBackendDataFrame works", {
    be <- backendInitialize(ChromBackendDataFrame(), chromData = test_df)

    new_ints <- lapply(test_df$rtime, function(z) z / 2)
    rtime(be) <- new_ints
    expect_identical(rtime(be), NumericList(new_ints, compress = FALSE))

    expect_error(rtime(be) <- 3, "has to be a list")
    expect_error(rtime(be) <- list(3, 2), "match the length")
    expect_error(rtime(be) <- list(3, 2, 4), "number of data points")
    new_ints[[2]] <- c(3.2, 1.2)
    expect_error(rtime(be) <- new_ints, "sorted increasingly")
})

test_that("selectChromVariables,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    res <- selectChromVariables(be, c("dataStorage", "msLevel"))

    df <- DataFrame(msLevel = 1:2, mz = c(2.3, 1.2),
                    other_col = 2)
    be <- backendInitialize(ChromBackendDataFrame(), df)

    res <- selectChromVariables(be, c("dataStorage", "other_col"))
    expect_equal(colnames(res@chromData), c("dataStorage", "other_col"))
    expect_equal(msLevel(res), c(NA_integer_, NA_integer_))

    res <- selectChromVariables(be, c("dataStorage", "mz"))
    expect_equal(colnames(res@chromData), c("dataStorage", "mz"))

    expect_error(selectChromVariables(be, "rtime"), "dataStorage is/are missing")
    expect_error(selectChromVariables(be, "something"),
                 "something not available")
})

test_that("$,$<-,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_identical(be$msLevel, integer())

    df <- DataFrame(msLevel = 1:2, mz = c(2.3, 1.2),
                    other_col = 2)
    be <- backendInitialize(be, df)
    expect_identical(be$msLevel, 1:2)
    expect_identical(be$other_col, c(2, 2))

    be$other_col <- 4
    expect_equal(be$other_col, c(4, 4))

    df$rtime <- list(1:3, 1:4)
    df$intensity <- list(c(3, 3, 3), c(4, 4, 4, 4))
    be <- backendInitialize(ChromBackendDataFrame(), df)
    be$intensity <- list(c(5, 5, 5), 1:4)
    expect_equal(be$intensity, NumericList(c(5, 5, 5), 1:4, compress = FALSE))
})

test_that("[,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_error(be[1])
    df <- DataFrame(scanIndex = 1:2, a = "a", b = "b")
    be <- backendInitialize(be, df)
    res <- be[1]
    expect_equal(be@chromData[1, ], res@chromData[1, ])
    res <- be[2]
    expect_equal(be@chromData[2, ], res@chromData[1, ])
    res <- be[2:1]
    expect_equal(be@chromData[2:1, ], res@chromData)

    res <- be[c(FALSE, FALSE)]
    expect_true(length(res) == 0)
    res <- be[c(FALSE, TRUE)]
    expect_equal(be@chromData[2, ], res@chromData[1, ])

    expect_error(be[TRUE], "match the length of")
    expect_error(be["a"], "does not have names")

    df <- DataFrame(scanIndex = c(1L, 2L, 1L, 2L),
                    file = c("a", "a", "b", "b"))
    be <- backendInitialize(be, df)
    dataStorage(be) <- c("1", "1", "2", "2")
    res <- be[3]
    expect_equal(dataStorage(res), "2")
    expect_equal(res@chromData$file, "b")

    res <- be[c(3, 1)]
    expect_equal(dataStorage(res), c("2", "1"))
    expect_equal(res@chromData$file, c("b", "a"))
})

test_that("filterDataOrigin,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterDataOrigin(be))
    df <- DataFrame(mz = as.numeric(1:8),
                    msLevel = 1L)
    be <- backendInitialize(ChromBackendDataFrame(), df)
    dataStorage(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    expect_true(length(filterDataOrigin(be, c(3, 1))) == 0)
    expect_equal(be, filterDataOrigin(be, NA_character_))

    dataOrigin(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterDataOrigin(be, c(3, 1))
    expect_equal(length(res), 6)
    expect_equal(dataOrigin(res), c("3", "3", "3", "1", "1", "1"))
    expect_equal(mz(res), c(6, 7, 8, 1, 2, 3))
    expect_equal(unique(res$dataStorage), c("3", "1"))

    res <- filterDataOrigin(be, c("2", "1"))
    expect_equal(length(res), 5)
    expect_equal(dataOrigin(res), c("2", "2", "1", "1", "1"))
    expect_equal(mz(res), c(4, 5, 1, 2, 3))

    res <- filterDataOrigin(be, 2)
    expect_equal(mz(res), c(4, 5))
    expect_equal(unique(res$dataStorage), "2")

    res <- filterDataOrigin(be, c(2, 3))
    expect_equal(mz(res), c(4, 5, 6, 7, 8))
    expect_equal(unique(res$dataStorage), c("2", "3"))

    res <- filterDataOrigin(be)
    expect_equal(res, be)

    df$dataOrigin <- "1"
    be <- backendInitialize(ChromBackendDataFrame(), df)
    res <- filterDataOrigin(be, 1L)
    expect_equal(res, be)
})

test_that("filterDataStorage,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterDataStorage(be))
    df <- DataFrame(mz = as.numeric(1:8),
                    msLevel = 1L)
    be <- backendInitialize(ChromBackendDataFrame(), df)
    dataStorage(be) <- c("1", "1", "1", "2", "2", "3", "3", "3")
    res <- filterDataStorage(be, c(3, 1))
    expect_equal(length(res), 6)
    expect_equal(dataStorage(res), c("3", "3", "3", "1", "1", "1"))
    expect_equal(mz(res), c(6, 7, 8, 1, 2, 3))
    expect_equal(unique(res$dataStorage), c("3", "1"))

    res <- filterDataStorage(be, c("2", "1"))
    expect_equal(length(res), 5)
    expect_equal(dataStorage(res), c("2", "2", "1", "1", "1"))
    expect_equal(mz(res), c(4, 5, 1, 2, 3))

    res <- filterDataStorage(be, 2)
    expect_equal(mz(res), c(4, 5))
    expect_equal(unique(res$dataStorage), "2")

    res <- filterDataStorage(be, c(2, 3))
    expect_equal(mz(res), c(4, 5, 6, 7, 8))
    expect_equal(unique(res$dataStorage), c("2", "3"))

    res <- filterDataStorage(be)
    expect_equal(res, be)
})

test_that("filterMsLevel,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterMsLevel(be))

    df <- DataFrame(msLevel = c(1L, 2L, 1L, 2L, 3L, 2L),
                    mz = as.numeric(1:6))
    be <- backendInitialize(ChromBackendDataFrame(), df)
    res <- filterMsLevel(be, 2L)
    expect_true(all(msLevel(res) == 2))
    expect_equal(mz(res), c(2, 4, 6))

    res <- filterMsLevel(be, c(3L, 2L))
    expect_equal(mz(res), c(2, 4, 5, 6))

    res <- filterMsLevel(be, c(3L, 5L))
    expect_equal(mz(res), 5)
})

test_that("filterMz,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterMz(be))

    df <- DataFrame(msLevel = 1:6,
                    mz = c(1, 2, 1, 1, 3, 4))
    be <- backendInitialize(ChromBackendDataFrame(), df)
    expect_equal(be, filterMz(be))
    res <- filterMz(be, 5)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    res <- filterMz(be, 1)
    expect_equal(msLevel(res), c(1L, 3L, 4L))
    expect_true(all(mz(res) == 1))

    res <- filterMz(be, c(2, 3))
    expect_equal(msLevel(res), c(2L, 5L))

    res <- filterMz(be, c(4, 10))
    expect_equal(msLevel(res), 6L)
})

test_that("filterPrecursorMz,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterPrecursorMz(be))

    df <- DataFrame(msLevel = 1:6,
                    precursorMzMin = c(1, 2, 1, 1, 3, 4),
                    precursorMzMax = c(1.1, 2.1, 1.1, 1.1, 3.1, 4.1))
    be <- backendInitialize(ChromBackendDataFrame(), df)
    expect_equal(be, filterPrecursorMz(be))
    res <- filterPrecursorMz(be, 5)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    res <- filterPrecursorMz(be, 1)
    expect_equal(msLevel(res), c(1L, 3L, 4L))
    expect_true(all(precursorMzMin(res) == 1))

    res <- filterPrecursorMz(be, c(2, 3))
    expect_equal(msLevel(res), c(2L, 5L))

    res <- filterPrecursorMz(be, c(4, 10))
    expect_equal(msLevel(res), 6L)

    ## precursorMzMin being NA.
    df <- DataFrame(msLevel = 1:6)
    be <- backendInitialize(be, df)
    res <- filterPrecursorMz(be, 1)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    be$precursorMzMin <- c(1, 2, 1, 1, 3, 4)
    res <- filterPrecursorMz(be, 1)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = 1:3,
                    precursorMz = c(1.1, 2.2, 3.3))
    be <- backendInitialize(be, df)
    res <- filterPrecursorMz(be, c(1, 2))
    expect_equal(msLevel(res), 1L)

    res <- filterPrecursorMz(be, c(2, 4))
    expect_equal(msLevel(res), 2:3)
})

test_that("filterProductMz,ChromBackendDataFrame works", {
    be <- ChromBackendDataFrame()
    expect_equal(be, filterProductMz(be))

    df <- DataFrame(msLevel = 1:6,
                    productMzMin = c(1, 2, 1, 1, 3, 4),
                    productMzMax = c(1.1, 2.1, 1.1, 1.1, 3.1, 4.1))
    be <- backendInitialize(ChromBackendDataFrame(), df)
    expect_equal(be, filterProductMz(be))
    res <- filterProductMz(be, 5)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    res <- filterProductMz(be, 1)
    expect_equal(msLevel(res), c(1L, 3L, 4L))
    expect_true(all(productMzMin(res) == 1))

    res <- filterProductMz(be, c(2, 3))
    expect_equal(msLevel(res), c(2L, 5L))

    res <- filterProductMz(be, c(4, 10))
    expect_equal(msLevel(res), 6L)

    ## productMzMin being NA.
    df <- DataFrame(msLevel = 1:6)
    be <- backendInitialize(be, df)
    res <- filterProductMz(be, 1)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    be$precursorMzMin <- c(1, 2, 1, 1, 3, 4)
    res <- filterProductMz(be, 1)
    expect_true(is(res, "ChromBackendDataFrame"))
    expect_true(length(res) == 0)

    df <- DataFrame(msLevel = 1:3,
                    productMz = c(1.1, 2.2, 3.3))
    be <- backendInitialize(be, df)
    res <- filterProductMz(be, c(1, 2))
    expect_equal(msLevel(res), 1L)

    res <- filterProductMz(be, c(2, 4))
    expect_equal(msLevel(res), 2:3)
})

test_that("split,ChromBackendDataFrame works", {
    chb <- mrm_mzr
    chbl <- split(chb, f = chb$dataStorage)
    expect_true(is(chbl[[1]], "ChromBackendDataFrame"))
    expect_identical(chb, chbl[[1]])
    expect_true(is(chbl[[1]]@chromData$polarity, "integer"))

    chbl <- split(chb, f = 1:length(chb))
    expect_true(is(chbl[[1]], "ChromBackendDataFrame"))
    expect_identical(chbl[[1]]$intensity[[1]], chb$intensity[[1]])

    chb2 <- backendMerge(chbl)
    expect_identical(chb2, chb)
})
