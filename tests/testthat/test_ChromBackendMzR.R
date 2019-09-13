test_that("initializeBackend,ChromBackendMzR works", {
    fl <- msdata::proteomics(full.names = TRUE)[1]
    expect_error(backendInitialize(ChromBackendMzR()), "Parameter 'files'")
    expect_error(backendInitialize(ChromBackendMzR(), files = 4),
                 "expected to be a character")
    expect_warning(
        expect_error(backendInitialize(ChromBackendMzR(), files = "some"),
                     "not found"))
    be <- backendInitialize(ChromBackendMzR(), files = fl)
    expect_true(validObject(be))
    expect_true(is(be, "ChromBackendMzR"))
    expect_equal(unique(be$dataStorage), normalizePath(fl))
    expect_equal(length(be), 138)
    expect_equal(be$chromIndex, 1:138)
    expect_equal(be@chromData$dataStorage, Rle(rep(normalizePath(fl), 138)))
    expect_true(isReadOnly(be))

    show(be)
    show(ChromBackendMzR())
})

test_that("chromData,chromData<-,ChromBackendMzR works", {
    be <- ChromBackendMzR()
    res <- chromData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)
    expect_true(all(names(Chromatograms:::.CHROMATOGRAMS_DATA_COLUMNS) %in%
                    colnames(res)))

    be <- mrm_mzr
    res <- .chrom_data_mzR(be, columns = c("msLevel", "chromIndex"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == length(be))
    expect_true(all(is.na(res$msLevel)))
    expect_true(is.integer(res$msLevel))
    expect_identical(res$chromIndex, 1:138)

    res <- chromData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == length(be))
    expect_true(all(is.na(res$msLevel)))
    expect_true(is.integer(res$msLevel))
    expect_true(is(res$intensity, "NumericList"))
    expect_true(is(res$rtime, "NumericList"))
    expect_true(all(names(Chromatograms:::.CHROMATOGRAMS_DATA_COLUMNS) %in%
                    colnames(res)))

    be$new_col <- "b"
    res <- chromData(be)
    any(colnames(res) == "new_col")
    expect_true(all(res$new_col == "b"))

    res$msLevel <- 2L
    expect_warning(chromData(be) <- res)
    expect_true(all(msLevel(be) == 2))
})

test_that("intensity,intensity<-,ChromBackendMzR works", {
    be <- ChromBackendMzR()
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == 0)
    expect_error(intensity(be) <- list(), "does not support replacing")

    be <- mrm_mzr
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
})

test_that("rtime,rtime<-,ChromaBackendMzR works", {
    be <- ChromBackendMzR()
    res <- rtime(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == 0)
    expect_error(rtime(be) <- list(), "does not support replacing")

    be <- mrm_mzr
    res <- rtime(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    be <- be[c(5, 2, 3)]
    res_2 <- rtime(be)
    expect_equal(res_2, res[c(5, 2, 3)])
})

test_that("$,ChromBackendMzR works", {
    be <- ChromBackendMzR()
    be$msLevel <- integer()
    expect_error(be$rtime <- list(), "not support replacing")
    expect_error(be$intensity <- list(), "not support replacing")

    be <- mrm_mzr
    be$msLevel <- 1L
    expect_true(all(msLevel(be) == 1))
    expect_error(msLevel(be) <- c(2L, 4L), "of length 138")
    be$msLevel <- 1:138
    expect_identical(be$msLevel, 1:138)
})
