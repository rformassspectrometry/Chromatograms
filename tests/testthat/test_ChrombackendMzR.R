## test backendInitialize

test_that("ChromBackendMzR works", {
    expect_true(isReadOnly(be_mzr))
    expect_false(be_mzr@inMemory)
    expect_false(identical(peaksData(be_mzr), be_mzr@peaksData))
    expect_true(identical(length(peaksData(be_mzr)), length((be_mzr@peaksData))))
    expect_true(identical(chromData(be_mzr),
                          fillCoreChromVariables(be_mzr@chromData)))
})

test_that("replacement method works", {
    tmp <- be_mzr
    peaksData(tmp)[[1]] <- peaksData(be_mzr)[[1]] + 1
    expect_false(identical(peaksData(tmp)[[1]], peaksData(be_mzr)[[1]]))
    expect_true(identical(peaksData(tmp)[[1]], peaksData(be_mzr)[[1]] + 1))
    expect_true(tmp@inMemory)
    expect_false(be_mzr@inMemory)

    cd <- chromData(tmp)
    cd$mz <-  1
    expect_true(!identical(chromData(tmp), cd))
    chromData(tmp) <- cd
    expect_equal(chromData(tmp), cd)
})

test_that("backendParallelFactor works", {
    expect_equal(levels(backendParallelFactor(be_mzr)), unique(dataOrigin(be_mzr)))
})
