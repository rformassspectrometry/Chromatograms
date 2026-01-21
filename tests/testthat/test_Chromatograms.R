test_that("Chromatograms works", {
    ## empty object
    expect_true(is(.backend(c_empty), "ChromBackendMemory"))
    expect_true(is(c_empty, "Chromatograms"))
    expect_true(processingChunkSize(c_empty)== Inf)
    expect_true(c_empty@version == "0.1")
    expect_identical(.processingQueue(c_empty), list())

    ## object with backend
    expect_true(is(.backend(c_full), "ChromBackendMemory"))
    expect_true(is(c_full, "Chromatograms"))
    expect_true(processingChunkSize(c_full) == Inf)
    expect_true(c_full@version == "0.1")
    expect_identical(.processingQueue(c_full), list())

    expect_equal(processingChunkSize(c_full), Inf)
    c_chunk <- c_full
    processingChunkSize(c_chunk) <- 2
    expect_equal(processingChunkSize(c_chunk), 2)
    expect_equal(levels(processingChunkFactor(c_chunk)), c("1", "2"))

    ## method with Spectra works
    c_sp <- Chromatograms(s[1:2])
    expect_true(is(.backend(c_sp), "ChromBackendSpectra"))
    expect_true(is(c_sp, "Chromatograms"))
    expect_true(processingChunkSize(c_sp) == Inf)
    expect_true(c_sp@version == "0.1")
    expect_identical(.processingQueue(c_sp), list())
})

test_that("Chromatograms constructor from Spectra works with all parameters", {
    ## Basic construction with defaults
    chr <- Chromatograms(s)
    expect_s4_class(chr, "Chromatograms")
    expect_s4_class(.backend(chr), "ChromBackendSpectra")
    expect_equal(length(chr), 3L)
    
    ## With summarize.method = "sum" (default)
    chr_sum <- Chromatograms(s, summarize.method = "sum")
    expect_s4_class(chr_sum, "Chromatograms")
    expect_identical(.backend(chr_sum)@summaryFun, sumi)
    
    ## With summarize.method = "max"
    chr_max <- Chromatograms(s, summarize.method = "max")
    expect_s4_class(chr_max, "Chromatograms")
    expect_identical(.backend(chr_max)@summaryFun, maxi)
    
    ## With empty chromData (should create default)
    chr_empty_cd <- Chromatograms(s, chromData = data.frame())
    expect_s4_class(chr_empty_cd, "Chromatograms")
    expect_true(nrow(chromData(chr_empty_cd)) > 0)
    expect_true(all(coreChromVariables() %in% colnames(chromData(chr_empty_cd))))
    
    ## With custom chromData
    custom_cd <- data.frame(
        msLevel = 1L,
        dataOrigin = unique(s$dataOrigin),
        customCol = "test"
    )
    chr_custom <- Chromatograms(s, chromData = custom_cd)
    expect_s4_class(chr_custom, "Chromatograms")
    expect_true("customCol" %in% colnames(chromData(chr_custom)))
    expect_equal(chromData(chr_custom)$customCol, "test")
    
    ## With custom factorize.by
    chr_factby <- Chromatograms(s, factorize.by = "dataOrigin")
    expect_s4_class(chr_factby, "Chromatograms")
    expect_true(all(chromData(chr_factby)$dataOrigin == 
                    chromData(chr_factby)$chromSpectraIndex))
    
    ## With spectraVariables
    chr_specvars <- Chromatograms(s, spectraVariables = c("polarity"))
    expect_s4_class(chr_specvars, "Chromatograms")
    if ("polarity" %in% Spectra::spectraVariables(s)) {
        expect_true("polarity" %in% colnames(chromData(chr_specvars)))
    }
})

test_that("Chromatograms constructor from ChromBackend works", {
    ## From ChromBackendMemory
    chr_mem <- Chromatograms(be)
    expect_s4_class(chr_mem, "Chromatograms")
    expect_s4_class(.backend(chr_mem), "ChromBackendMemory")
    expect_equal(length(chr_mem), length(be))
    
    ## From ChromBackendMzR
    chr_mzr <- Chromatograms(be_mzr)
    expect_s4_class(chr_mzr, "Chromatograms")
    expect_s4_class(.backend(chr_mzr), "ChromBackendMzR")
    expect_equal(length(chr_mzr), length(be_mzr))
    
    ## From ChromBackendSpectra
    chr_spec <- Chromatograms(be_sp)
    expect_s4_class(chr_spec, "Chromatograms")
    expect_s4_class(.backend(chr_spec), "ChromBackendSpectra")
    expect_equal(length(chr_spec), length(be_sp))
    
    ## With processingQueue
    pq <- list(function(x) x)
    chr_pq <- Chromatograms(be, processingQueue = pq)
    expect_equal(length(.processingQueue(chr_pq)), 1)
})

test_that("Chromatograms constructor handles edge cases", {
    ## Empty Spectra
    empty_s <- Spectra()
    chr_empty <- Chromatograms(empty_s)
    expect_s4_class(chr_empty, "Chromatograms")
    expect_equal(length(chr_empty), 0)
    
    ## Missing object (creates empty ChromBackendMemory)
    chr_missing <- Chromatograms()
    expect_s4_class(chr_missing, "Chromatograms")
    expect_s4_class(.backend(chr_missing), "ChromBackendMemory")
    expect_equal(length(chr_missing), 0)
})

test_that("show, Chromatograms - ChromBackendMemory works", {
    expect_output(show(c_full), "ChromBackendMemory")
    res <- c_full
    res@processing <- c("a", "b", "c", "d")
    expect_output(show(res), "1 more processings")
    res@processingQueue <- list("a", "b", "c", "d")
    expect_output(show(res), "4 processing step")
})

test_that("show, Chromatograms - ChromBackendMzR works", {
    expect_output(show(c_mzr), "ChromBackendMzR")
    res <- c_mzr
    res@processing <- c("a", "b", "c", "d")
    expect_output(show(res), "1 more processings")
    res@processingQueue <- list("a", "b", "c", "d")
    expect_output(show(res), "4 processing step")
})


test_that("setBackend works correctly", {
    c_mzr_new <- setBackend(c_mzr, backend = ChromBackendMemory())
    expect_s4_class(.backend(c_mzr_new), "ChromBackendMemory")
    expect_identical(chromData(c_mzr_new), chromData(c_mzr))
    expect_identical(peaksData(c_mzr_new), peaksData(c_mzr))
    expect_identical(peaksData(c_mzr_new), peaksData(c_mzr))

    processingChunkSize(c_mzr) <- 100
    f <- processingChunkFactor(c_mzr)
    expect_true(length(levels(f)) > 1)
    c_mzr_new <- setBackend(c_mzr, backend = ChromBackendMemory(), f = f)
    expect_s4_class(.backend(c_mzr_new), "ChromBackendMemory")
    expect_identical(chromData(c_mzr_new), chromData(c_mzr))
    expect_identical(peaksData(c_mzr_new), peaksData(c_mzr))
    expect_identical(.peaksData(c_mzr_new), peaksData(c_mzr))

    expect_error(
        setBackend(c_mzr, backend = ChromBackendMzR()),
        "does not support"
    )

    c_sp_new <- setBackend(c_sp, backend = ChromBackendMemory())
    expect_true(!all(c("rtMin", "rtMax") %in% colnames(chromData(c_sp_new))))
})

test_that("$ works correctly", {
    expect_identical(msLevel(c_full), c_full$msLevel)
    expect_identical(chromIndex(c_mzr), c_mzr$chromIndex)
    expect_identical(intensity(c_full), c_full$intensity)
    expect_identical(intensity(c_mzr), c_mzr$intensity)
    tmp <- c_full
    tmp$msLevel <- c(2L, 2L, 3L)
    expect_identical(msLevel(tmp), c(2L, 2L, 3L))
    tmp$intensity <- lapply(tmp$intensity, function(x) x + 10)
    expect_false(identical(intensity(tmp), intensity(c_full)))
})

test_that("[ works correctly", {
    c_sub <- c_full[1:2]
    expect_true(is(c_sub, "Chromatograms"))
    expect_equal(nrow(chromData(c_sub)), 2)
    expect_equal(length(peaksData(c_sub)), 2)
    expect_error(c_full[1:2, 1], "by columns is not")

    c_sub <- c_full[1]
    expect_equal(c_sub, c_sub[])
})

test_that("[[ works properly", {
    expect_error(c_full[[1]], "character")
    expect_error(c_full[["test", 1]], "not supported")
    expect_error(c_full[["test"]], "No variable")
    expect_equal(c_full[["msLevel"]], msLevel(c_full))
    expect_equal(c_full[["msLevel"]], .backend(c_full)[["msLevel"]])

    ## replace
    expect_error(c_full[[1]] <- 1, "character defining the chromatogram")
    expect_error(c_full[["msLevel", 4]] <- 1, "not supported")
    expect_error(c_full[["test"]] <- 1, "No variable")
    repet <- c_full
    repet[["msLevel"]] <- rep(2L, length(repet))
    expect_false(identical(msLevel(repet), msLevel(c_full)))
})

test_that("factorize() works", {
    tmp <- c_sp
    tmp$msLevel <- c(1L, 2L, 3L)
    idx_before <- chromSpectraIndex(.backend(tmp))
    tmp <- factorize(tmp)
    idx_after <- chromSpectraIndex(.backend(tmp))
    expect_false(identical(idx_before, idx_after))
})

test_that("chromExtract, Chromatograms works correctly", {
    peak.table <- data.frame(
        msLevel = 1L,
        dataOrigin = "mem1",
        rtMin = 12.5,
        rtMax = 14.0
    )
    res <- chromExtract(c_full, peak.table, by = c("msLevel", "dataOrigin"))
    expect_s4_class(res, "Chromatograms")
    expect_true(length(res) == nrow(peak.table))

    pk_nomatch <- transform(peak.table, dataOrigin = "no_such_sample")
    expect_error(chromExtract(c_full, pk_nomatch,
                              by = c("msLevel", "dataOrigin")),
                 "Some combinations in")

})

