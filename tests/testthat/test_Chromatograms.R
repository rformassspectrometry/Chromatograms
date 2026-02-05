test_that("Chromatograms works", {
  ## empty object
  expect_true(is(.backend(c_empty), "ChromBackendMemory"))
  expect_true(is(c_empty, "Chromatograms"))
  expect_true(processingChunkSize(c_empty) == Inf)
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
  expect_true(all(
    names(coreChromVariables()) %in%
      chromVariables(chr_empty_cd)
  ))

  ## With custom chromData
  custom_cd <- data.frame(
    msLevel = 1L,
    dataOrigin = unique(s$dataOrigin),
    customCol = "test"
  )
  chr_custom <- Chromatograms(s, chromData = custom_cd)
  expect_s4_class(chr_custom, "Chromatograms")
  expect_true("customCol" %in% colnames(chromData(chr_custom)))
  expect_equal(
    chromData(chr_custom)$customCol,
    rep("test", length(unique(s$dataOrigin)))
  )

  ## With custom factorize.by
  chr_factby <- Chromatograms(s, factorize.by = "dataOrigin")
  expect_s4_class(chr_factby, "Chromatograms")
  expect_true(all(
    chromData(chr_factby)$dataOrigin == chromData(chr_factby)$chromSpectraIndex
  ))

  ## With spectraVariables
  expect_error(
    Chromatograms(s, spectraVariables = "polarity"),
    "must already exist in the chromData"
  )

  chr_specvars <- Chromatograms(s, spectraVariables = c("precursorMz"))
  expect_s4_class(chr_specvars, "Chromatograms")
  if ("precursorMz" %in% Spectra::spectraVariables(s)) {
    expect_true("precursorMz" %in% colnames(chromData(chr_specvars)))
  }

  ## spectraVariables should replace all-NA columns in provided chromData
  sp_small <- Spectra::Spectra(S4Vectors::DataFrame(
    rtime = c(1, 2),
    msLevel = c(1L, 1L),
    dataOrigin = c("A", "A"),
    polarity = c(1L, 1L)
  ))
  cd_na <- data.frame(msLevel = 1L, dataOrigin = "A", polarity = NA_integer_)
  chr_specvars_replace <- Chromatograms(
    sp_small,
    chromData = cd_na,
    spectraVariables = c("polarity")
  )
  expect_s4_class(chr_specvars_replace, "Chromatograms")
  expect_identical(chromData(chr_specvars_replace)$polarity, 1L)
})

test_that("Chromatograms constructor from ChromBackend works", {
  ## From ChromBackendMemory - already initialized
  chr_mem <- Chromatograms(be)
  expect_s4_class(chr_mem, "Chromatograms")
  expect_s4_class(.backend(chr_mem), "ChromBackendMemory")
  expect_equal(length(chr_mem), length(be))

  ## From empty ChromBackendMemory with chromData and peaksData parameters
  cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem1", "mem1")
  )
  pdata <- list(
    data.frame(
      rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
      intensity = c(100, 250, 400, 300, 150)
    ),
    data.frame(rtime = numeric(), intensity = numeric()),
    data.frame(
      rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
      intensity = c(80, 500, 1200, 600, 120)
    )
  )
  chr <- Chromatograms(
    ChromBackendMemory(),
    chromData = cdata,
    peaksData = pdata
  )
  expect_s4_class(chr, "Chromatograms")
  expect_s4_class(.backend(chr), "ChromBackendMemory")
  expect_equal(length(chr), 3L)
  expect_identical(chromData(chr)$mz, cdata$mz)
  expect_identical(peaksData(chr), pdata)

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
  pq <- list(ProcessingStep(
    "smooth",
    list(method = "SavitzkyGolay", halfWindowSize = 2L)
  ))
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

test_that("setBackend preserves peaksData order with duplicate chromSpectraIndex", {
  ## This test reproduces a bug where setBackend would scramble peaksData order
  ## when multiple chromatograms share the same chromSpectraIndex (e.g., EICs
  ## with different mz windows from the same file). The bug occurred because
  ## peaksData() was returning data in factor level order instead of row order.
  sp <- Spectra(DataFrame(
    mz = NumericList(
      c(100, 101, 102), c(100, 101, 102), c(100, 101, 102),
      c(100, 101, 102), c(100, 101, 102), c(100, 101, 102),
      compress = FALSE
    ),
    intensity = NumericList(
      c(10, 20, 30), c(15, 25, 35), c(20, 30, 40),
      c(100, 200, 300), c(150, 250, 350), c(200, 300, 400),
      compress = FALSE
    ),
    rtime = c(1, 2, 3, 4, 5, 6),
    msLevel = rep(1L, 6),
    dataOrigin = c("A", "A", "A", "B", "B", "B")
  ))
  
  ## Create multiple EICs from same spectra (sharing chromSpectraIndex)
  custom_chromData <- data.frame(
    msLevel = rep(1L, 4),
    mzMin = c(100, 100, 100, 100),
    mzMax = c(101, 102, 101, 102),
    rtMin = c(1, 1, 4, 4),
    rtMax = c(3, 3, 6, 6),
    dataOrigin = c("A", "A", "B", "B")
  )
  
  chr <- Chromatograms(sp, chromData = custom_chromData)
  
  ## Get peaksData before setBackend
  pd_before <- peaksData(chr)
  
  ## Convert to memory backend
  chr_mem <- setBackend(chr, ChromBackendMemory())
  
  ## Get peaksData after setBackend
  pd_after <- peaksData(chr_mem)
  
  ## peaksData should be identical before and after setBackend
  expect_identical(length(pd_before), length(pd_after))
  
  ## Each element should have the same intensity values
  for (i in seq_along(pd_before)) {
    expect_identical(pd_before[[i]]$intensity, pd_after[[i]]$intensity)
    expect_identical(pd_before[[i]]$rtime, pd_after[[i]]$rtime)
  }
  
  ## Verify subsetting after setBackend works correctly
  ## This tests the scenario where subsetting after setBackend returns wrong data
  chr_sub <- chr_mem[2]  # Get second chromatogram (1_A with mzMax=102)
  pd_sub <- peaksData(chr_sub)
  expect_identical(pd_sub[[1]]$intensity, pd_before[[2]]$intensity)
  
  chr_sub4 <- chr_mem[4]  # Get fourth chromatogram (1_B with mzMax=102)
  pd_sub4 <- peaksData(chr_sub4)
  expect_identical(pd_sub4[[1]]$intensity, pd_before[[4]]$intensity)
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
  expect_error(
    chromExtract(c_full, pk_nomatch, by = c("msLevel", "dataOrigin")),
    "Some combinations in"
  )
})
