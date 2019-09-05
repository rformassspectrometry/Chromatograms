test_that(".valid_chrom_backend_files_exist", {
    expect_match(.valid_chrom_backend_files_exist(c("a", "b")), "a, b not")
    tmpf <- tempfile()
    write("hello", file = tmpf)
    expect_null(.valid_chrom_backend_files_exist(tmpf))
    expect_null(.valid_chrom_backend_files_exist(character()))
    expect_null(.valid_chrom_backend_files_exist(NA_character_))
})

test_that(".mzR_chrom_header works", {
    fl <- msdata::proteomics(full.names = TRUE)[1]
    expect_error(.mzR_chrom_header(), "should have length 1")
    hdr <- .mzR_chrom_header(fl)
    expect_true(inherits(hdr, "DataFrame"))
    expect_equal(hdr$chromIndex, 1:nrow(hdr))
    ## first chromatogram is TIC
    expect_true(is.na(hdr$precursorMz[1]))

    fl <- dir(system.file("cdf", package = "msdata"), full.names = TRUE)
    expect_warning(hdr <- .mzR_chrom_header(fl))
    expect_true(nrow(hdr) == 0)
    expect_true(inherits(hdr, "DataFrame"))
})

test_that(".mzR_chromatograms works", {
    fl <- msdata::proteomics(full.names = TRUE)[1]
    expect_error(.mzR_chromatograms(), "should have length 1")
    res <- .mzR_chromatograms(fl, 2)
    expect_true(is.list(res))
    expect_equal(colnames(res[[1]]), c("rtime", "intensity"))
    expect_true(length(res) == 1)

    res <- .mzR_chromatograms(fl, 1:10)
    expect_true(is.list(res))
    expect_equal(colnames(res[[1]]), c("rtime", "intensity"))
    expect_true(length(res) == 10)
    expect_true(is.matrix(res[[1]]))
    expect_true(is.numeric(res[[1]][, "rtime"]))
    expect_true(is.numeric(res[[1]][, "intensity"]))
})

test_that(".chrom_data_mzR works", {
    ## Tested in test_ChromBackendMzR.R chromData,chromData<-
})

test_that(".rtime_intensity_pairs_mzR works", {
    fl <- msdata::proteomics(full.names = TRUE)[1]
    be <- ChromBackendMzR()
    res <- as.list(be)
    expect_true(is(res, "list"))
    expect_true(length(res) == 0)

    be <- mrm_mzr
    res <- .rtime_intensity_pairs_mzR(be)
    expect_true(is(res, "list"))
    expect_true(is.matrix(res[[2]]))
    expect_true(all(colnames(res[[2]]) == c("rtime", "intensity")))
})
