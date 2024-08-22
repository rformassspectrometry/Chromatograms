library(S4Vectors)

test_that("ChromBackend methods throw errors", {
    setClass("DummyBackend",
             contains = "ChromBackend")
    dm <- new("DummyBackend")

    expect_error(dm[1], "Not implemented for ")
    expect_error(dm$a, "Not implemented for ")
    expect_error(dm$a <- "a", "Not implemented for ")
    expect_error(backendMerge(dm), "Not implemented for ")
    expect_error(chromData(dm), "Not implemented for ")
    expect_error(chromData(dm) <- DataFrame(), "Not implemented for ")
    expect_error(peaksData(dm), "Not implemented for ")
    expect_error(peaksData(dm) <- list(), "Not implemented for ")
    expect_error(selectChromVariables(dm), "Not implemented for ")

    expect_true(isReadOnly(dm))
    expect_equal(backendParallelFactor(dm), factor())
})


