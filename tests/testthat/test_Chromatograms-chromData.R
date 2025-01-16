test_that("Chromatograms, chromData, chromData<- works", {
    expect_equal(chromData(c_full), fillCoreChromVariables(c_full@backend@chromData))
    cdata2 <- cdata
    cdata2$mz <- cdata2$mz + 1
    res <- c_full
    chromData(res) <- cdata2
    expect_equal(chromData(res), fillCoreChromVariables(cdata2))
    expect_true(all(chromData(res)$mz != fillCoreChromVariables(cdata)$mz))
})

test_that("Chromatograms, chromIndex, chromIndex<- works", {
    expect_equal(chromIndex(c_full), c(1, 2, 3))
    res <- c_full
    chromIndex(res) <- c(3L, 3L, 1L)
    expect_equal(chromIndex(res), c(3, 3, 1))
    expect_true(all(chromIndex(res) != c(1, 2, 3)))
})

test_that("Chromatograms, collisionEnergy, collisionEnergy<- works", {
    res <- c_full
    chromData(res)$collisionEnergy <- c(35.0, 40.0, 45.0)
    expect_equal(collisionEnergy(res), c(35.0, 40.0, 45.0))
    collisionEnergy(res) <- c(50.0, 55.0, 60.0)
    expect_equal(collisionEnergy(res), c(50.0, 55.0, 60.0))
})

test_that("Chromatograms, dataOrigin, dataOrigin<- works", {
    res <- c_full
    chromData(res)$dataOrigin <- c("file1", "file2", "file3")
    expect_equal(dataOrigin(res), c("file1", "file2", "file3"))
    dataOrigin(res) <- c("source1", "source2", "source3")
    expect_equal(dataOrigin(res), c("source1", "source2", "source3"))
})

test_that("Chromatograms, dataStorage, dataStorage<- works", {
    res <- c_full
    chromData(res)$dataStorage <- c("disk", "memory", "disk")
    expect_equal(dataStorage(res), c("disk", "memory", "disk"))
    dataStorage(res) <- c("memory", "memory", "disk")
    expect_equal(dataStorage(res), c("memory", "memory", "disk"))
})

test_that("Chromatograms, msLevel, msLevel<- works", {
    res <- c_full
    expect_equal(msLevel(res), c(1L, 1L, 1L))
    msLevel(res) <- c(2L, 3L, 1L)
    expect_equal(msLevel(res), c(2L, 3L, 1L))
})

test_that("Chromatograms, mz, mz<- works", {
    res <- c_full
    expect_equal(mz(res), c(112.2, 123.3, 134.4))
    mz(res) <- c(115.0, 125.0, 135.0)
    expect_equal(mz(res), c(115.0, 125.0, 135.0))
})

test_that("Chromatograms, mzMin, mzMin<- works", {
    res <- c_full
    chromData(res)$mzMin <- c(110.0, 120.0, 130.0)
    expect_equal(mzMin(res), c(110.0, 120.0, 130.0))
    mzMin(res) <- c(111.0, 121.0, 131.0)
    expect_equal(mzMin(res), c(111.0, 121.0, 131.0))
})

test_that("Chromatograms, mzMax, mzMax<- works", {
    res <- c_full
    chromData(res)$mzMax <- c(114.0, 124.0, 134.0)
    expect_equal(mzMax(res), c(114.0, 124.0, 134.0))
    mzMax(res) <- c(116.0, 126.0, 136.0)
    expect_equal(mzMax(res), c(116.0, 126.0, 136.0))
})

test_that("Chromatograms, length works", {
    expect_equal(length(c_full), 3)
    expect_equal(length(c_empty), 0)
})

test_that("Chromatograms, precursorMz, precursorMz<- works", {
    res <- c_full
    chromData(res)$precursorMz <- c(400.0, 500.0, 600.0)
    expect_equal(precursorMz(res), c(400.0, 500.0, 600.0))
    precursorMz(res) <- c(450.0, 550.0, 650.0)
    expect_equal(precursorMz(res), c(450.0, 550.0, 650.0))
})

test_that("Chromatograms, precursorMzMin, precursorMzMin<- works", {
    res <- c_full
    chromData(res)$precursorMzMin <- c(395.0, 495.0, 595.0)
    expect_equal(precursorMzMin(res), c(395.0, 495.0, 595.0))
    precursorMzMin(res) <- c(390.0, 490.0, 590.0)
    expect_equal(precursorMzMin(res), c(390.0, 490.0, 590.0))
})

test_that("Chromatograms, precursorMzMax, precursorMzMax<- works", {
    res <- c_full
    chromData(res)$precursorMzMax <- c(405.0, 505.0, 605.0)
    expect_equal(precursorMzMax(res), c(405.0, 505.0, 605.0))
    precursorMzMax(res) <- c(410.0, 510.0, 610.0)
    expect_equal(precursorMzMax(res), c(410.0, 510.0, 610.0))
})

test_that("Chromatograms, productMz, productMz<- works", {
    res <- c_full
    chromData(res)$productMz <- c(100.0, 150.0, 200.0)
    expect_equal(productMz(res), c(100.0, 150.0, 200.0))
    productMz(res) <- c(110.0, 160.0, 210.0)
    expect_equal(productMz(res), c(110.0, 160.0, 210.0))
})

test_that("Chromatograms, productMzMin, productMzMin<- works", {
    res <- c_full
    chromData(res)$productMzMin <- c(95.0, 145.0, 195.0)
    expect_equal(productMzMin(res), c(95.0, 145.0, 195.0))
    productMzMin(res) <- c(90.0, 140.0, 190.0)
    expect_equal(productMzMin(res), c(90.0, 140.0, 190.0))
})

test_that("Chromatograms, productMzMax, productMzMax<- works", {
    res <- c_full
    chromData(res)$productMzMax <- c(105.0, 155.0, 205.0)
    expect_equal(productMzMax(res), c(105.0, 155.0, 205.0))
    productMzMax(res) <- c(110.0, 160.0, 210.0)
    expect_equal(productMzMax(res), c(110.0, 160.0, 210.0))
})


test_that("filterChromData handles various edge cases", {
    expect_identical(
        filterChromData(be_cd, variables = c("mz"), ranges = numeric(),
                        match = "any"),
        be_cd
    )
    res <- filterChromData(be_cd, variables = c("mz"), ranges = c(100, 200),
                           match = "all")
    expect_identical(res, be_cd)

    res <- filterChromData(be_cd, variables = c("mz"), ranges = c(500, 600),
                           match = "any")
    expect_equal(nrow(chromData(res)), 0)
    expect_error(
        filterChromData(be_cd, variables = c("mz", "chromIndex"),
                        ranges = c(100, 200), match = "any"),
        "be twice the length of the "
    )
    expect_error(
        filterChromData(be_cd, variables = c("mz"),
                        ranges = c("a", "b"), match = "any"),
        "filterChromData only support filtering for numerical"
    )
    expect_error(
        filterChromData(be_cd, variables = c("nonExistentVar"),
                        ranges = c(100, 200), match = "any"),
        " not available"
    )
    res <- filterChromData(c_empty, variables = c("mz"),
                           ranges = c(100, 200), match = "any")
    expect_equal(nrow(chromData(res)), 0)
    res <- filterChromData(be_cd, variables = c("mz", "chromIndex"),
                           ranges = c(134, 150, 1, 2), match = "any")
    expect_true(nrow(chromData(res)) > 0)
    res <- filterChromData(be_cd, variables = c("mz", "chromIndex"),
                           ranges = c(134, 150, 1, 2), match = "all")
    expect_true(nrow(chromData(res)) <= nrow(chromData(be_cd)))
    res <- filterChromData(be_cd, variables = c("mz", "chromIndex"),
                           ranges = c(100, 200, 1, 3), match = "any")
    expect_equal(nrow(chromData(res)), 3)
    res <- filterChromData(be_cd, variables = c("mz", "chromIndex"),
                           ranges = c(500, 600, 10, 20), match = "all")
    expect_equal(nrow(chromData(res)), 0)
    res <- filterChromData(be_cd, variables = c("mz"), ranges = c(134, 150),
                           match = "any", keep = FALSE)
    expect_equal(nrow(chromData(res)), 2)
    res <- filterChromData(be_cd, variables = c("mz"), ranges = c(120, 130), match = "any")
    expect_equal(nrow(chromData(res)), 1)
})


