test_that("ChromBackend methods throw errors", {
    setClass("DummyBackend",
             contains = "ChromBackend")
    dm <- new("DummyBackend")

    expect_error(dm[1], "Not implemented for ")
    expect_error(dm$a, "Not implemented for ")
    expect_error(dm$a <- "a", "Not implemented for ")
    expect_error(backendMerge(dm), "Not implemented for ")
    expect_error(chromData(dm), "Not implemented for ")
    expect_error(chromData(dm) <- data.frame(), "Not implemented for ")
    expect_error(peaksData(dm), "Not implemented for ")
    expect_error(peaksData(dm) <- list(), "Not implemented for ")
    expect_true(!isReadOnly(dm))
    expect_equal(backendParallelFactor(dm), factor())

    expect_true(is(reset(dm), class(dm)))

    # accessor function with default but dependent on chromData() and
    # therefore not implemented
    expect_true(is(backendInitialize(dm), class(dm)))
    expect_error(chromIndex(dm), "Not implemented for ")
    expect_error(chromIndex(dm) <- 1, "Not implemented for ")
    expect_error(collisionEnergy(dm), "Not implemented for ")
    expect_error(collisionEnergy(dm) <- 1, "Not implemented for ")
    expect_error(dataOrigin(dm), "Not implemented for ")
    expect_error(dataOrigin(dm) <- "a", "Not implemented for ")
    expect_error(dataStorage(dm), "Not implemented for ")
    expect_error(dataStorage(dm) <- "a", "Not implemented for ")
    expect_error(msLevel(dm), "Not implemented for ")
    expect_error(msLevel(dm) <- 1, "Not implemented for ")
    expect_error(mz(dm), "Not implemented for ")
    expect_error(mz(dm) <- 1, "Not implemented for ")
    expect_error(mzMax(dm), "Not implemented for ")
    expect_error(mzMax(dm) <- 1, "Not implemented for ")
    expect_error(mzMin(dm), "Not implemented for ")
    expect_error(mzMin(dm) <- 1, "Not implemented for ")
    expect_error(precursorMz(dm), "Not implemented for ")
    expect_error(precursorMz(dm) <- 1, "Not implemented for ")
    expect_error(precursorMzMax(dm), "Not implemented for ")
    expect_error(precursorMzMax(dm) <- 1, "Not implemented for ")
    expect_error(precursorMzMin(dm), "Not implemented for ")
    expect_error(precursorMzMin(dm) <- 1, "Not implemented for ")
    expect_error(productMz(dm), "Not implemented for ")
    expect_error(productMz(dm) <- 1, "Not implemented for ")
    expect_error(productMzMax(dm), "Not implemented for ")
    expect_error(productMzMax(dm) <- 1, "Not implemented for ")
    expect_error(productMzMin(dm), "Not implemented for ")
    expect_error(productMzMin(dm) <- 1, "Not implemented for ")
    expect_error(rtime(dm), "Not implemented for ")
    expect_error(rtime(dm) <- c(1,1,1), "Not implemented for ")

    expect_identical(chromVariables(dm), names(coreChromVariables()))
    expect_identical(peaksVariables(dm), names(corePeaksVariables()))

    expect_error(dm[[123]], "is supposed to be a character")
    expect_error(dm[["character", "character2"]], "is not supported")
    expect_error(dm[[123]] <- 123, "is supposed to be a character")
    expect_error(dm[["character", "character2"]] <- 123, "is not supported")
    expect_identical(filterPeaksData(dm, ranges = c()), dm)
    expect_identical(filterPeaksData(dm, ranges = c(1,2), variables = c()), dm)
    expect_error(filterPeaksData(dm, ranges = c("a","3"), variables = c("a")),
                 "only support")
    expect_error(filterPeaksData(dm, ranges = c(1,2), variables = c("a")),
                 "One or more")
    expect_error(filterPeaksData(dm, ranges = c(1,2), variables = c(33)),
                 "needs to be of type")
    expect_error(filterPeaksData(dm, ranges = c(1, 2),
                                 variables = c("intensity", "rtime")),
                 "needs to be twice")

})


