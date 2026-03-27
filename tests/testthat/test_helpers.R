test_that(".df_combine works as expected", {
  combined_backend <- .df_combine(list(be, be_cd))
  expect_equal(nrow(combined_backend@chromData), nrow(cdata) + nrow(cdata))
  expect_equal(
    length(.peaksData(combined_backend)),
    length(.peaksData(be)) + length(.peaksData(be_cd))
  )

  expect_equal(.df_combine(list(be)), be)
  incompatible_data <- list(
    data.frame(
      rtime = c(10.0, 12.0),
      intensity = c(200, 150),
      other_col = c("test", "test")
    ),
    data.frame(
      rtime = c(30.1, 31.2),
      intensity = c(110, 90),
      other_col = c("test", "test")
    )
  )
  be_incompatible <- backendInitialize(
    be_empty,
    chromData = cdata,
    peaksData = incompatible_data
  )

  expect_error(
    .df_combine(c(be, be_incompatible)),
    "Provided objects have different sets of peak variables."
  )

  setClass("DummyBackend", contains = "ChromBackend")
  dm <- new("DummyBackend")

  expect_error(
    .df_combine(c(be, dm)),
    "merge backends of the same type"
  )
})


test_that(".filter_ranges helper function works correctly", {
  query <- data.frame(
    mz = c(112.2, 123.3, 134.4),
    chromIndex = c(1L, 2L, 3L)
  )
  ranges <- numeric()
  expect_error(.filter_ranges(query, ranges, match = "any"), "needs to be")
  ranges <- c(1, 2, 3)
  expect_error(
    .filter_ranges(query, ranges, match = "any"),
    "Length of 'ranges'"
  )

  # All ranges match (full data retained)
  ranges <- c(100, 200, 1, 3)
  res <- .filter_ranges(query, ranges, match = "all")
  expect_equal(res, c(1, 2, 3))
  res <- .filter_ranges(query, ranges, match = "any")
  expect_equal(res, c(1, 2, 3))

  # No matches (empty result)
  ranges <- c(500, 600, 4, 5)
  res <- .filter_ranges(query, ranges, match = "any")
  expect_equal(res, integer(0))

  # Partially overlapping ranges
  ranges <- c(120, 130, 4, 5)
  res <- .filter_ranges(query, ranges, match = "any")
  expect_equal(res, 2)
  res <- .filter_ranges(query, ranges, match = "all")
  expect_equal(res, integer(0))

  # Single row edge case
  query_single <- data.frame(mz = 125.0, chromIndex = 2L)
  ranges <- c(120, 130, 1, 3)
  res <- .filter_ranges(query_single, ranges, match = "all")
  expect_equal(res, 1)
  res <- .filter_ranges(query_single, ranges, match = "any")
  expect_equal(res, 1)

  # Edge case: empty query
  query_empty <- data.frame(mz = numeric(), chromIndex = integer())
  ranges <- c(100, 200, 1, 3)
  res <- .filter_ranges(query_empty, ranges, match = "all")
  expect_equal(res, integer(0))
  res <- .filter_ranges(query_empty, ranges, match = "any")
  expect_equal(res, integer(0))
})

test_that(".check_column_order_and_types works", {
  df_valid <- data.frame(
    rtime = c(1, 2, 3),
    intensity = c(10, 20, 30)
  )
  df_invalid_order <- data.frame(
    intensity = c(10, 20, 30),
    rtime = c(1, 2, 3)
  )
  df_invalid_type <- data.frame(
    rtime = c("a", "b", "c"),
    intensity = c(10, 20, 30)
  )
  df_empty <- data.frame(rtime = numeric(), intensity = numeric())

  expect_null(.check_column_order_and_types(
    df_valid,
    names(.CORE_PEAKS_VARIABLES),
    .CORE_PEAKS_VARIABLES
  ))
  expect_equal(
    .check_column_order_and_types(
      df_invalid_order,
      names(.CORE_PEAKS_VARIABLES),
      .CORE_PEAKS_VARIABLES
    ),
    "Columns should be in the order 'rtime', 'intensity'."
  )
  expect_equal(
    .check_column_order_and_types(
      df_invalid_type,
      names(.CORE_PEAKS_VARIABLES),
      .CORE_PEAKS_VARIABLES
    ),
    "The peaksData variable(s) rtime have the wrong data type."
  )
  expect_null(.check_column_order_and_types(
    df_empty,
    names(.CORE_PEAKS_VARIABLES),
    .CORE_PEAKS_VARIABLES
  ))
})


test_that(".check_rtime works", {
  df_valid <- data.frame(
    rtime = c(1, 2, 3),
    intensity = c(10, 20, 30)
  )
  df_invalid_na <- data.frame(
    rtime = c(1, NA, 3),
    intensity = c(10, 20, 30)
  )
  df_invalid_increasing <- data.frame(
    rtime = c(3, 2, 1),
    intensity = c(10, 20, 30)
  )
  df_empty <- data.frame(rtime = numeric(), intensity = numeric())

  expect_null(.check_rtime(df_valid))
  expect_equal(
    .check_rtime(df_invalid_na),
    "'rtime' column contains NA values."
  )
  expect_equal(
    .check_rtime(df_invalid_increasing),
    "'rtime' column is not strictly increasing."
  )
  expect_null(.check_rtime(df_empty))
})

test_that(".validate_entry works", {
  df_valid <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
  df_empty <- data.frame(rtime = numeric(), intensity = numeric())

  expect_null(.validate_entry(
    df_valid,
    1,
    names(.CORE_PEAKS_VARIABLES),
    .CORE_PEAKS_VARIABLES
  ))
  expect_null(.validate_entry(
    df_empty,
    3,
    names(.CORE_PEAKS_VARIABLES),
    .CORE_PEAKS_VARIABLES
  ))
  expect_equal(
    .validate_entry(
      c(1, 2, 3),
      1,
      names(.CORE_PEAKS_VARIABLES),
      .CORE_PEAKS_VARIABLES
    ),
    "Entry 1: all 'peaksData' entries should be of class 'data.frame'"
  )
})

test_that(".run_process_queue ChromBackendMemory work", {
  result <- .run_process_queue(
    c_empty@backend,
    f = processingChunkFactor(c_empty),
    queue = c_empty@processingQueue
  )
  expect_equal(result, c_empty@backend)

  result <- .run_process_queue(
    c_full@backend,
    f = processingChunkFactor(c_full),
    queue = c_full@processingQueue
  )
  expect_equal(result, c_full@backend)

  c_queued <- filterPeaksData(
    c_full,
    variables = c("rtime"),
    ranges = c(12.5, 45.5)
  )
  c_queued <- filterPeaksData(
    c_queued,
    variables = c("intensity"),
    ranges = c(100, 200)
  )

  # this test for f = factor() and queue >1
  result <- .run_process_queue(
    c_queued@backend,
    f = processingChunkFactor(c_queued),
    queue = c_queued@processingQueue
  )
  expect_true(inherits(result, "ChromBackend"))
  peaks_result <- peaksData(result)
  expect_equal(length(peaks_result), length(peaksData(c_full@backend)))
  expect_false(identical(peaks_result, peaksData(c_full@backend)))

  f <- factor(c(1, 1, 2))
  c_queued <- filterPeaksData(
    c_full,
    variables = c("rtime"),
    ranges = c(12.5, 45.5)
  )
  c_queued <- filterPeaksData(
    c_queued,
    variables = c("intensity"),
    ranges = c(100, 200)
  )

  result <- .run_process_queue(
    c_queued@backend,
    f = f,
    queue = c_queued@processingQueue
  )

  expect_true(inherits(result, "ChromBackend"))
  expect_equal(length(peaksData(result)), length(peaksData(c_full@backend)))

  split_data <- split(c_full@backend, f)
  expect_equal(length(split_data), length(levels(f)))

  f <- factor(c(1, 2))
  tmp <- c_full
  tmp@processingQueue <- list(1)
  expect_error(
    .run_process_queue(tmp@backend, f = f, queue = tmp@processingQueue),
    "length 'f' has to be equal to the length of 'object'"
  )

  f <- c(1, 2, 3)
  expect_error(
    .run_process_queue(tmp@backend, f = f, queue = tmp@processingQueue),
    "f must be a factor"
  )
})

test_that(".run_processing_queue, ChromBackendMzr work", {
  ## without factor and queue == 1
  c_queued <- filterPeaksData(
    c_mzr,
    variables = c("rtime"),
    ranges = c(12.5, 25.5),
    keep = FALSE
  )

  result1 <- .run_process_queue(
    c_queued@backend,
    f = processingChunkFactor(c_queued),
    queue = c_queued@processingQueue
  )

  expect_true(result1@inMemory)
  expect_false(c_queued@backend@inMemory)

  expect_true(inherits(result1, "ChromBackendMzR"))
  expect_false(identical(
    lengths(rtime(result1)),
    lengths(rtime(c_mzr@backend))
  ))
  expect_equal(length(peaksData(result1)), length(peaksData(c_mzr@backend)))

  ## with levels(factor) > 1 and queue == 1
  processingChunkSize(c_queued) <- 100
  f <- processingChunkFactor(c_queued) # > 1
  result2 <- .run_process_queue(
    c_queued@backend,
    f = f,
    queue = c_queued@processingQueue
  )
  expect_true(inherits(result2, "ChromBackendMzR"))
  expect_equal(length(peaksData(result2)), length(peaksData(c_mzr@backend)))
  expect_false(identical(
    lengths(rtime(result2)),
    lengths(rtime(c_mzr@backend))
  ))
  expect_true(result2@inMemory)
  expect_false(c_queued@backend@inMemory)
  expect_identical(result1, result2)

  ## without factor and queue > 1
  c_queued <- filterPeaksData(
    c_queued,
    variables = c("intensity"),
    ranges = c(45, 50)
  )
  result3 <- .run_process_queue(
    c_queued@backend,
    f = factor(),
    queue = c_queued@processingQueue
  )

  expect_true(inherits(result3, "ChromBackendMzR"))
  expect_equal(length(peaksData(result3)), length(peaksData(c_mzr@backend)))
  expect_false(identical(
    lengths(rtime(result3)),
    lengths(rtime(c_mzr@backend))
  ))
  expect_true(result3@inMemory)
  expect_false(c_queued@backend@inMemory)

  ## with factor and queue > 1
  f <- processingChunkFactor(c_queued)
  result4 <- .run_process_queue(
    c_queued@backend,
    f = f,
    queue = c_queued@processingQueue
  )
  expect_true(inherits(result4, "ChromBackendMzR"))
  expect_equal(
    length(peaksData(result4)),
    length(peaksData(c_mzr@backend))
  )
  expect_false(identical(
    lengths(rtime(result4)),
    lengths(rtime(c_mzr@backend))
  ))
  expect_true(result4@inMemory)
  expect_false(c_queued@backend@inMemory)
  expect_identical(result3, result4)
})

test_that(".valid_processing_queue works correctly", {
  expect_null(.valid_processing_queue(list()))
  valid_queue <- list(new("ProcessingStep"))
  expect_null(.valid_processing_queue(valid_queue))

  invalid_queue <- list("not_a_processing_step")
  expect_error(
    .valid_processing_queue(invalid_queue),
    "'processingQueue' should only contain ProcessingStep objects."
  )
})

test_that("ensure_rt_mz_columns correctly handles mz and rt columns", {
  spectra <- s
  spectra_f <- factor(
    do.call(
      paste,
      c(
        as.list(Spectra::spectraData(s)[, c("msLevel", "dataOrigin")]),
        sep = "_"
      )
    )
  )
  levs <- levels(spectra_f)
  chrom_data <- data.frame(msLevel = c(1, 2, 3), chromSpectraIndex = levs[1:3])
  chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
  expect_equal(chrom_data$mzMin, c(-Inf, -Inf, -Inf))
  expect_equal(chrom_data$mzMax, c(Inf, Inf, Inf))

  chrom_data <- data.frame(mzMin = c(100), chromSpectraIndex = levs[1])
  expect_error(
    .ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
    "must be present if one is provided."
  )

  chrom_data <- data.frame(mzMax = c(200), chromSpectraIndex = levs[1])
  expect_error(
    .ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
    "must be present if one is provided."
  )
  chrom_data <- data.frame(msLevel = c(1, 2, 3), chromSpectraIndex = levs[1:3])
  chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
  s_plit <- split(spectra, spectra_f)
  expect_equal(chrom_data$rtMin[[1]], min(s_plit[[1]]$rtime, na.rm = TRUE))
  expect_equal(chrom_data$rtMax[[1]], max(s_plit[[1]]$rtime, na.rm = TRUE))

  chrom_data <- data.frame(rtMin = c(10), chromSpectraIndex = levs[1])
  expect_error(
    .ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
    " must be present if one is provided."
  )
  chrom_data <- data.frame(rtMax = c(50), chromSpectraIndex = levs[1])
  expect_error(
    .ensure_rt_mz_columns(chrom_data, spectra, spectra_f),
    "must be present if one is provided."
  )

  chrom_data <- data.frame(
    mzMin = c(100),
    mzMax = c(200),
    rtMin = c(10),
    rtMax = c(50),
    chromSpectraIndex = levs[1]
  )
  chrom_data <- .ensure_rt_mz_columns(chrom_data, spectra, spectra_f)
  expect_equal(chrom_data$mzMin, 100)
  expect_equal(chrom_data$mzMax, 200)
  expect_equal(chrom_data$rtMin, 10)
  expect_equal(chrom_data$rtMax, 50)
})

test_that(".validate_chromExtract_input works correctly", {
  cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    dataOrigin = c("A", "B", "C")
  )
  be <- backendInitialize(new("ChromBackendMemory"), chromData = cdata)

  peak_tbl <- data.frame(
    rtMin = c(1, 2, 3),
    rtMax = c(5, 6, 7),
    msLevel = c(1L, 1L, 1L),
    dataOrigin = c("A", "B", "C")
  )

  # should pass
  expect_silent(
    .validate_chromExtract_input(be, peak_tbl, by = c("msLevel", "dataOrigin"))
  )

  # missing required column
  bad_tbl <- peak_tbl[, !names(peak_tbl) %in% "rtMax", drop = FALSE]
  expect_error(
    .validate_chromExtract_input(be, bad_tbl, by = c("msLevel", "dataOrigin")),
    "must contain columns"
  )

  # NA in rtMin
  bad_tbl2 <- peak_tbl
  bad_tbl2$rtMin[1] <- NA
  expect_error(
    .validate_chromExtract_input(be, bad_tbl2, by = c("msLevel", "dataOrigin")),
    "cannot contain NA"
  )

  # missing 'by' columns in chromData
  bad_be <- backendInitialize(
    new("ChromBackendMemory"),
    chromData = cdata[, "msLevel", drop = FALSE]
  )
  expect_error(
    .validate_chromExtract_input(
      bad_be,
      peak_tbl,
      by = c("msLevel", "dataOrigin")
    ),
    "must be present"
  )

  ## unique
  expect_error(
    .validate_chromExtract_input(bad_be, peak_tbl, by = "msLevel"),
    "must uniquely identify rows"
  )
})

test_that(".match_chromdata_peaktable aligns correctly", {
  tmp_cdata <- data.frame(
    msLevel = c(1L, 1L, 2L),
    dataOrigin = c("A", "B", "A")
  )
  tmp <- backendInitialize(new("ChromBackendMemory"), chromData = tmp_cdata)

  peak_tbl <- data.frame(
    msLevel = c(1L, 2L),
    dataOrigin = c("A", "A"),
    rtMin = c(1, 2),
    rtMax = c(5, 6)
  )

  matched <- .match_chromdata_peaktable(
    tmp,
    peak_tbl,
    by = c("msLevel", "dataOrigin")
  )

  # Expect a subset of object
  expect_s4_class(matched$object, "ChromBackendMemory")
  expect_equal(length(matched$chrom_keys), nrow(.chromData(matched$object)))

  # Check factor levels alignment
  expect_true(all(levels(matched$peak_keys) %in% levels(matched$chrom_keys)))

  # missing key should error
  bad_tbl <- data.frame(
    msLevel = 3L,
    dataOrigin = "Z",
    rtMin = 1,
    rtMax = 2
  )
  expect_error(
    .match_chromdata_peaktable(tmp, bad_tbl, by = c("msLevel", "dataOrigin")),
    "do not exist"
  )
})

test_that(".check_overl_columns warns correctly", {
  tmp_cdata <- data.frame(
    msLevel = 1L,
    dataOrigin = "X",
    mz = 100,
    extracol = "info"
  )
  tmp <- backendInitialize(new("ChromBackendMemory"), chromData = tmp_cdata)

  peak_tbl <- data.frame(
    rtMin = 1,
    rtMax = 2,
    mzMin = 99,
    mzMax = 101,
    msLevel = 1L,
    dataOrigin = "X",
    mz = 123,
    extracol = "test"
  )

  req_cols <- c("rtMin", "rtMax", "mzMin", "mzMax", "msLevel", "dataOrigin")

  expect_warning(
    overl <- .check_overl_columns(tmp, peak_tbl, req_cols),
    "already exist"
  )

  # overlapping should include "mz"
  expect_true(all(c("mz", "extracol") %in% names(peak_tbl)[overl]))
})

test_that(".impute() works correctly and without warnings", {
  # Base signal with gaps
  x <- c(1:5, NA, 7:10, NA, 12:15, rep(NA, 2), 18:20)

  ## linear
  expect_silent({
    res_lin <- .impute(x, method = "linear")
  })
  expect_false(anyNA(res_lin))
  expect_true(all(diff(res_lin) > 0)) # still increasing

  ## spline
  expect_silent({
    res_spl <- .impute(x, method = "spline")
  })
  expect_false(anyNA(res_spl))
  expect_equal(length(res_spl), length(x))

  ## Gaussian
  expect_silent(
    res_gauss <- .impute(x, method = "gaussian", window = 2, sd = 1)
  )
  expect_false(anyNA(res_gauss))
  expect_equal(length(res_gauss), length(x))

  ## loess
  expect_warning({
    res_loess <- .impute(x, method = "loess", span = 0.3)
    "could not fill all NAs"
  })
  expect_false(anyNA(res_loess))
  expect_equal(length(res_loess), length(x))

  ## Consecutive NAs
  x_na <- c(1, 2, NA, NA, 5, 6, 7, 8, NA, NA, 11, 12)
  expect_silent({
    res_consec <- .impute(x_na, method = "linear")
  })
  expect_false(anyNA(res_consec))
  expect_equal(length(res_consec), length(x_na))

  ## No NA
  x_nomiss <- 1:10
  expect_silent({
    res_nomiss <- .impute(x_nomiss, method = "spline")
  })
  expect_identical(res_nomiss, x_nomiss)

  ## all Nas returns NA
  x_allna <- rep(NA_real_, 8)
  expect_silent({
    res_allna <- .impute(x_allna, method = "gaussian")
  })
  expect_true(all(is.na(res_allna)))
})

test_that(".map_spectra_vars() correctly maps spectra variables", {
  vars <- c("scanIndex", "mtbls_id")

  mapped <- .map_spectra_vars(be_sp, vars)

  expect_s4_class(mapped, "ChromBackendSpectra")
  cd <- chromData(mapped)

  expect_true(all(vars %in% names(cd)))

  expect_true(all(vapply(
    cd$mtbls_id,
    function(x) length(unique(x)) == 1,
    logical(1)
  )))

  expect_true(all(vapply(cd$scanIndex, function(x) length(x) >= 1, logical(1))))

  expect_identical(
    rownames(cd),
    rownames(chromData(be_sp))
  )
  expect_error(
    .map_spectra_vars(be_sp, c("scanIndex", "fake_var")),
    "must exist in the Spectra"
  )

  cd_names <- names(chromData(be_sp))
  fake_var <- cd_names[1L] # pick a real chromData column
  expect_error(
    .map_spectra_vars(be_sp, fake_var),
    "must already exist in the chromData"
  )
})

test_that(".prepare_spectra_input extracts and classifies correctly", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), c(100, 200, 300),
                     compress = FALSE),
    intensity = NumericList(c(10, 20, 30), c(15, 25, 35), c(20, 30, 40),
                            compress = FALSE),
    rtime = c(1, 2, 3),
    msLevel = rep(1L, 3),
    dataOrigin = rep("A", 3)
  ))

  ## EIC case (finite mz range)
  cd <- data.frame(rtMin = 1, rtMax = 3, mzMin = 100, mzMax = 200)
  prep <- .prepare_spectra_input(cd, sp, sumi)
  expect_equal(prep$n_cd, 1L)
  expect_equal(prep$n_valid, 3L)
  expect_true(prep$same_rt)
  expect_false(prep$all_tic)
  expect_true(prep$use_cumsum)
  expect_equal(prep$valid_rt, c(1, 2, 3))
  expect_length(prep$valid_mz, 3L)
  expect_length(prep$valid_int, 3L)
  ## Each mz/int entry should be a plain numeric vector
  expect_type(prep$valid_mz[[1]], "double")
  expect_type(prep$valid_int[[1]], "double")

  ## TIC case (infinite mz ranges)
  cd_tic <- data.frame(rtMin = 1, rtMax = 3, mzMin = -Inf, mzMax = Inf)
  prep_tic <- .prepare_spectra_input(cd_tic, sp, sumi)
  expect_true(prep_tic$all_tic)

  ## maxi function → use_cumsum = FALSE
  prep_max <- .prepare_spectra_input(cd, sp, maxi)
  expect_false(prep_max$use_cumsum)

  ## Multiple chromatograms with same rt
  cd_multi <- data.frame(
    rtMin = c(1, 1), rtMax = c(3, 3),
    mzMin = c(100, 200), mzMax = c(200, 300)
  )
  prep_multi <- .prepare_spectra_input(cd_multi, sp, sumi)
  expect_equal(prep_multi$n_cd, 2L)
  expect_true(prep_multi$same_rt)
  expect_false(prep_multi$all_tic)

  ## Different rt ranges → same_rt = FALSE
  cd_diff <- data.frame(
    rtMin = c(1, 2), rtMax = c(2, 3),
    mzMin = c(100, 100), mzMax = c(200, 200)
  )
  prep_diff <- .prepare_spectra_input(cd_diff, sp, sumi)
  expect_false(prep_diff$same_rt)

  ## Narrow rt filter reduces n_valid
  cd_narrow <- data.frame(rtMin = 2, rtMax = 2, mzMin = 100, mzMax = 200)
  prep_narrow <- .prepare_spectra_input(cd_narrow, sp, sumi)
  expect_equal(prep_narrow$n_valid, 1L)
  expect_equal(prep_narrow$valid_rt, 2)
})

test_that(".prepare_spectra_input handles duplicated retention times", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), c(100, 200), c(100, 200),
                     compress = FALSE),
    intensity = NumericList(c(10, 20), c(15, 25), c(30, 40), c(35, 45),
                            compress = FALSE),
    rtime = c(1, 1, 2, 3),
    msLevel = rep(1L, 4),
    dataOrigin = rep("A", 4)
  ))
  cd <- data.frame(rtMin = 1, rtMax = 3, mzMin = -Inf, mzMax = Inf)
  prep <- .prepare_spectra_input(cd, sp, sumi)
  ## Only unique rtimes should be kept
  expect_equal(prep$n_valid, 3L)
  expect_equal(prep$valid_rt, c(1, 2, 3))
})


test_that(".spectrum_intensities TIC case returns replicated fun(int_j)", {
  mz_j <- c(100, 200, 300)
  int_j <- c(10, 20, 30)
  ## sumi: sum = 60, replicated for 3 chromatograms
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(-Inf, -Inf, -Inf), mzMaxs = c(Inf, Inf, Inf),
    fun = sumi, all_tic = TRUE, use_cumsum = TRUE
  )
  expect_equal(res, rep(60, 3))
  ## maxi: max = 30
  res_max <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(-Inf, -Inf), mzMaxs = c(Inf, Inf),
    fun = maxi, all_tic = TRUE, use_cumsum = FALSE
  )
  expect_equal(res_max, rep(30, 2))
})

test_that(".spectrum_intensities cumsum EIC case", {
  mz_j <- c(100, 200, 300)
  int_j <- c(10, 20, 30)
  ## Two ranges: mz 100-200 → 10+20=30, mz 200-300 → 20+30=50
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 200), mzMaxs = c(200, 300),
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_equal(res, c(30, 50))
})

test_that(".spectrum_intensities generic (maxi) case", {
  mz_j <- c(100, 200, 300)
  int_j <- c(10, 20, 30)
  ## mz 100-200 → max(10,20)=20, mz 200-300 → max(20,30)=30
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 200), mzMaxs = c(200, 300),
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_equal(res, c(20, 30))
})

test_that(".spectrum_intensities empty mz returns NA", {
  ## cumsum path
  res_cs <- .spectrum_intensities(numeric(0), numeric(0),
    mzMins = 100, mzMaxs = 200,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_true(is.na(res_cs))
  ## generic path
  res_gen <- .spectrum_intensities(numeric(0), numeric(0),
    mzMins = 100, mzMaxs = 200,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_true(is.na(res_gen))
})

test_that(".spectrum_intensities no mz in range returns NA", {
  mz_j <- c(500, 600)
  int_j <- c(10, 20)
  res_cs <- .spectrum_intensities(mz_j, int_j,
    mzMins = 100, mzMaxs = 200,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_true(is.na(res_cs))
  res_gen <- .spectrum_intensities(mz_j, int_j,
    mzMins = 100, mzMaxs = 200,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_true(is.na(res_gen))
})

test_that(".spectrum_intensities TIC with NAs uses fun directly", {
  mz_j <- c(100, 200, 300)
  int_j <- c(10, NA, 30)
  ## sumi (na.rm): 10 + 30 = 40 (MsCoreUtils::sumi removes NAs)
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = -Inf, mzMaxs = Inf,
    fun = sumi, all_tic = TRUE, use_cumsum = TRUE
  )
  expect_equal(res, sumi(int_j))
  ## maxi (na.rm): max = 30
  res_max <- .spectrum_intensities(mz_j, int_j,
    mzMins = -Inf, mzMaxs = Inf,
    fun = maxi, all_tic = TRUE, use_cumsum = FALSE
  )
  expect_equal(res_max, maxi(int_j))
})

test_that(".spectrum_intensities cumsum with NAs strips NA peaks", {
  mz_j <- c(100, 200, 300, 400)
  int_j <- c(10, NA, 30, NA)
  ## After stripping NAs: mz = c(100, 300), int = c(10, 30)
  ## mz 100-300: 10+30=40, mz 300-400: 30
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 300), mzMaxs = c(300, 400),
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_equal(res, c(40, 30))
})

test_that(".spectrum_intensities generic (maxi) with NAs", {
  mz_j <- c(100, 200, 300, 400)
  int_j <- c(10, NA, 30, NA)
  ## mz 100-300: maxi(10, NA, 30) = 30
  res <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 300), mzMaxs = c(300, 400),
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_equal(res[1], 30) # max(10, NA, 30) via maxi = 30
  expect_equal(res[2], 30) # max(30, NA) via maxi = 30
})

test_that(".spectrum_intensities all-NA spectrum returns NA", {
  mz_j <- c(100, 200, 300)
  int_j <- c(NA_real_, NA_real_, NA_real_)
  ## cumsum: .strip_na_peaks returns NULL → NA
  res_cs <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 200), mzMaxs = c(200, 300),
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_true(all(is.na(res_cs)))
  ## generic: maxi(NA, NA, NA) → NA
  res_gen <- .spectrum_intensities(mz_j, int_j,
    mzMins = c(100, 200), mzMaxs = c(200, 300),
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_true(all(is.na(res_gen)))
})


test_that(".build_intensity_matrix TIC/BPC fast path", {
  mz_list <- list(c(100, 200, 300), c(100, 200, 300))
  int_list <- list(c(10, 20, 30), c(15, 25, 35))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1:2,
    mzMins = c(-Inf, -Inf), mzMaxs = c(Inf, Inf), n_cd = 2L,
    fun = sumi, all_tic = TRUE, use_cumsum = TRUE
  )
  expect_equal(nrow(int_mat), 2L)
  expect_equal(ncol(int_mat), 2L)
  ## TIC: sum of all intensities per spectrum, same for all chroms
  expect_equal(int_mat[1, 1], 60)  # 10+20+30
  expect_equal(int_mat[1, 2], 60)
  expect_equal(int_mat[2, 1], 75)  # 15+25+35
  expect_equal(int_mat[2, 2], 75)
})

test_that(".build_intensity_matrix cumsum EIC case", {
  mz_list <- list(c(100, 200, 300), c(100, 200, 300))
  int_list <- list(c(10, 20, 30), c(15, 25, 35))
  ## 2 chromatograms: mz 100-200 and mz 200-300
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1:2,
    mzMins = c(100, 200), mzMaxs = c(200, 300), n_cd = 2L,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_equal(nrow(int_mat), 2L)
  expect_equal(ncol(int_mat), 2L)
  ## Chrom 1 (mz 100-200): spec1 = 10+20=30, spec2 = 15+25=40
  expect_equal(int_mat[1, 1], 30)
  expect_equal(int_mat[2, 1], 40)
  ## Chrom 2 (mz 200-300): spec1 = 20+30=50, spec2 = 25+35=60
  expect_equal(int_mat[1, 2], 50)
  expect_equal(int_mat[2, 2], 60)
})

test_that(".build_intensity_matrix generic function (maxi)", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(10, 20, 30))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = c(100, 200), mzMaxs = c(200, 300), n_cd = 2L,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_equal(int_mat[1, 1], 20)  # max(10, 20)
  expect_equal(int_mat[1, 2], 30)  # max(20, 30)
})

test_that(".build_intensity_matrix handles empty mz and no match", {
  mz_list <- list(numeric(0), c(100, 200), c(500, 600))
  int_list <- list(numeric(0), c(10, 20), c(50, 60))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1:3,
    mzMins = 100, mzMaxs = 200, n_cd = 1L,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_true(is.na(int_mat[1, 1]))   # empty mz → NA
  expect_equal(int_mat[2, 1], 30)     # 10+20
  expect_true(is.na(int_mat[3, 1]))   # no mz in range → NA
})


test_that(".build_intensity_matrix BPC multi-range (maxi)", {
  ## Multiple chromatograms with different mz ranges
  mz_list <- list(c(100, 150, 200, 250, 300),
                  c(100, 150, 200, 250, 300))
  int_list <- list(c(10, 50, 20, 80, 30),
                   c(15, 60, 25, 70, 35))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1:2,
    mzMins = c(100, 200), mzMaxs = c(200, 300), n_cd = 2L,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  ## Chrom 1 (mz 100-200): spec1 max(10,50,20)=50, spec2 max(15,60,25)=60
  expect_equal(int_mat[1, 1], 50)
  expect_equal(int_mat[2, 1], 60)
  ## Chrom 2 (mz 200-300): spec1 max(20,80,30)=80, spec2 max(25,70,35)=70
  expect_equal(int_mat[1, 2], 80)
  expect_equal(int_mat[2, 2], 70)
})

test_that(".build_intensity_matrix BPC single-element range", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(10, 20, 30))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = c(200), mzMaxs = c(200), n_cd = 1L,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_equal(int_mat[1, 1], 20)
})

test_that(".build_intensity_matrix BPC handles NA intensities", {
  mz_list <- list(c(100, 200, 300, 400))
  int_list <- list(c(10, NA, 30, NA))
  ## Sparse table path with NAs
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = c(100, 300), mzMaxs = c(300, 400), n_cd = 2L,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  ## mz 100-300: non-NA peaks are 10, 30 → max = 30

  expect_equal(int_mat[1, 1], 30)
  ## mz 300-400: non-NA peaks are 30 → max = 30
  expect_equal(int_mat[1, 2], 30)
})

test_that(".build_intensity_matrix cumsum handles NA intensities", {
  mz_list <- list(c(100, 200, 300, 400))
  int_list <- list(c(10, NA, 30, NA))
  int_mat <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = c(100, 300), mzMaxs = c(300, 400), n_cd = 2L,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  ## mz 100-300: non-NA are 10, 30 → sum = 40
  expect_equal(int_mat[1, 1], 40)
  ## mz 300-400: non-NA is 30 → sum = 30
  expect_equal(int_mat[1, 2], 30)
})

test_that(".build_intensity_matrix all-NA spectrum returns NA", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(NA_real_, NA_real_, NA_real_))
  ## cumsum path
  res_sum <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = 100, mzMaxs = 300, n_cd = 1L,
    fun = sumi, all_tic = FALSE, use_cumsum = TRUE
  )
  expect_true(is.na(res_sum[1, 1]))
  ## BPC path
  res_max <- .build_intensity_matrix(
    mz_list, int_list, kept = 1L,
    mzMins = 100, mzMaxs = 300, n_cd = 1L,
    fun = maxi, all_tic = FALSE, use_cumsum = FALSE
  )
  expect_true(is.na(res_max[1, 1]))
})

test_that(".compute_chrom_intensities TIC case", {
  mz_list <- list(c(100, 200), c(100, 200, 300))
  int_list <- list(c(10, 20), c(15, 25, 35))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1:2,
    mz_lo = -Inf, mz_hi = Inf, fun = sumi, use_cumsum = TRUE
  )
  expect_equal(res, c(30, 75))
})

test_that(".compute_chrom_intensities cumsum case", {
  mz_list <- list(c(100, 200, 300), c(100, 200, 300))
  int_list <- list(c(10, 20, 30), c(15, 25, 35))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1:2,
    mz_lo = 100, mz_hi = 200, fun = sumi, use_cumsum = TRUE
  )
  expect_equal(res, c(30, 40))
})

test_that(".compute_chrom_intensities generic function (maxi)", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(10, 20, 30))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1L,
    mz_lo = 100, mz_hi = 300, fun = maxi, use_cumsum = FALSE
  )
  expect_equal(res, 30)
})

test_that(".compute_chrom_intensities handles empty mz and no match", {
  mz_list <- list(numeric(0), c(500, 600))
  int_list <- list(numeric(0), c(10, 20))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1:2,
    mz_lo = 100, mz_hi = 200, fun = sumi, use_cumsum = TRUE
  )
  expect_true(is.na(res[1]))
  expect_true(is.na(res[2]))
})

test_that(".compute_chrom_intensities TIC with NAs", {
  mz_list <- list(c(100, 200), c(100, 200, 300))
  int_list <- list(c(10, NA), c(NA, 25, 35))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1:2,
    mz_lo = -Inf, mz_hi = Inf, fun = sumi, use_cumsum = TRUE
  )
  ## sumi removes NAs: spec1=10, spec2=60
  expect_equal(res, c(sumi(c(10, NA)), sumi(c(NA, 25, 35))))
})

test_that(".compute_chrom_intensities cumsum with NAs", {
  mz_list <- list(c(100, 200, 300), c(100, 200, 300))
  int_list <- list(c(10, NA, 30), c(NA, 25, 35))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1:2,
    mz_lo = 100, mz_hi = 200, fun = sumi, use_cumsum = TRUE
  )
  ## After stripping NAs: spec1 mz=c(100,300) int=c(10,30) → mz 100-200 = 10
  ## spec2 mz=c(200,300) int=c(25,35) → mz 100-200: only mz=200 → 25
  expect_equal(res[1], 10)
  expect_equal(res[2], 25)
})

test_that(".compute_chrom_intensities generic (maxi) with NAs", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(10, NA, 30))
  res <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1L,
    mz_lo = 100, mz_hi = 300, fun = maxi, use_cumsum = FALSE
  )
  ## maxi(10, NA, 30) = 30  (maxi ignores NAs)
  expect_equal(res, 30)
})

test_that(".compute_chrom_intensities all-NA spectrum returns NA", {
  mz_list <- list(c(100, 200, 300))
  int_list <- list(c(NA_real_, NA_real_, NA_real_))
  ## cumsum: .strip_na_peaks → NULL → NA
  res_sum <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1L,
    mz_lo = 100, mz_hi = 300, fun = sumi, use_cumsum = TRUE
  )
  expect_true(is.na(res_sum))
  ## generic path
  res_max <- .compute_chrom_intensities(
    mz_list, int_list, kept = 1L,
    mz_lo = 100, mz_hi = 300, fun = maxi, use_cumsum = FALSE
  )
  expect_true(is.na(res_max))
})


test_that(".process_peaks_data same-rt TIC case", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, 20, 30), c(15, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(
    rtMin = c(1, 1), rtMax = c(2, 2),
    mzMin = c(-Inf, -Inf), mzMax = c(Inf, Inf)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  expect_length(res, 2L)
  expect_true(all(vapply(res, is.data.frame, logical(1))))
  expect_equal(res[[1]]$rtime, c(1, 2))
  expect_equal(res[[1]]$intensity, c(60, 75))
  ## Both TIC with same rt → identical
  expect_equal(res[[1]], res[[2]])
})

test_that(".process_peaks_data same-rt EIC (mz filtering)", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, 20, 30), c(15, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(
    rtMin = c(1, 1), rtMax = c(2, 2),
    mzMin = c(100, 200), mzMax = c(200, 300)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  expect_equal(res[[1]]$intensity, c(30, 40))
  expect_equal(res[[2]]$intensity, c(50, 60))
})

test_that(".process_peaks_data different-rt case", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, 20), c(15, 25), c(30, 40), compress = FALSE),
    rtime = c(1, 2, 3),
    msLevel = rep(1L, 3),
    dataOrigin = rep("A", 3)
  ))
  cd <- data.frame(
    rtMin = c(1, 2), rtMax = c(2, 3),
    mzMin = c(-Inf, -Inf), mzMax = c(Inf, Inf)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  expect_equal(res[[1]]$rtime, c(1, 2))
  expect_equal(res[[1]]$intensity, c(30, 40))
  expect_equal(res[[2]]$rtime, c(2, 3))
  expect_equal(res[[2]]$intensity, c(40, 70))
})

test_that(".process_peaks_data with maxi function", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, 20, 30), c(15, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(rtMin = 1, rtMax = 2, mzMin = 100, mzMax = 300)
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), maxi, FALSE)
  expect_equal(res[[1]]$intensity, c(30, 35))
})

test_that(".process_peaks_data drop=TRUE returns vectors", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, 20), compress = FALSE),
    rtime = 1,
    msLevel = 1L,
    dataOrigin = "A"
  ))
  cd <- data.frame(rtMin = 1, rtMax = 1, mzMin = -Inf, mzMax = Inf)
  res <- .process_peaks_data(cd, sp, "intensity", sumi, drop = TRUE)
  expect_type(res[[1]], "double")
  expect_equal(res[[1]], 30)
})

test_that(".process_peaks_data handles empty rt range", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, 20), compress = FALSE),
    rtime = 1,
    msLevel = 1L,
    dataOrigin = "A"
  ))
  cd <- data.frame(rtMin = 5, rtMax = 10, mzMin = -Inf, mzMax = Inf)
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  expect_equal(nrow(res[[1]]), 0L)
})

test_that(".process_peaks_data rtime-only request", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, 20), c(15, 25), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(rtMin = 1, rtMax = 2, mzMin = 100, mzMax = 200)
  res <- .process_peaks_data(cd, sp, "rtime", sumi, drop = TRUE)
  expect_equal(res[[1]], c(1, 2))
})

test_that(".process_peaks_data returns correct results via ChromBackendSpectra", {
  ## Integration test using the actual fixture backend
  pd <- peaksData(be_sp)
  expect_true(is.list(pd))
  expect_true(all(vapply(pd, is.data.frame, logical(1))))
  expect_true(all(vapply(pd, function(d) all(c("rtime", "intensity") %in% names(d)),
                         logical(1))))
  ## rtime should be monotonically non-decreasing in each chromatogram
  for (i in seq_along(pd)) {
    if (nrow(pd[[i]]) > 1)
      expect_true(all(diff(pd[[i]]$rtime) >= 0))
  }
})

test_that(".process_peaks_data different-rt EIC with cumsum", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 150, 200), c(100, 150, 200), c(100, 150, 200),
                     compress = FALSE),
    intensity = NumericList(c(10, 5, 20), c(15, 8, 25), c(30, 12, 40),
                            compress = FALSE),
    rtime = c(1, 2, 3),
    msLevel = rep(1L, 3),
    dataOrigin = rep("A", 3)
  ))
  cd <- data.frame(
    rtMin = c(1, 2), rtMax = c(2, 3),
    mzMin = c(100, 100), mzMax = c(150, 200)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  ## Chrom 1: rt 1-2, mz 100-150 → spec1: 10+5=15, spec2: 15+8=23
  expect_equal(res[[1]]$intensity, c(15, 23))
  ## Chrom 2: rt 2-3, mz 100-200 → spec2: 15+8+25=48, spec3: 30+12+40=82
  expect_equal(res[[2]]$intensity, c(48, 82))
})

test_that(".process_peaks_data different-rt EIC with maxi", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, 20, 30), c(15, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(
    rtMin = c(1, 2), rtMax = c(1, 2),
    mzMin = c(100, 100), mzMax = c(200, 300)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), maxi, FALSE)
  ## Chrom 1: rt 1 only, mz 100-200 → max(10, 20) = 20
  expect_equal(res[[1]]$intensity, 20)
  ## Chrom 2: rt 2 only, mz 100-300 → max(15, 25, 35) = 35
  expect_equal(res[[2]]$intensity, 35)
})

test_that(".process_peaks_data same-rt TIC with NAs in intensities", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, NA, 30), c(NA, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(
    rtMin = c(1, 1), rtMax = c(2, 2),
    mzMin = c(-Inf, -Inf), mzMax = c(Inf, Inf)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  ## TIC: sumi ignores NAs → spec1: 10+30=40, spec2: 25+35=60
  expect_equal(res[[1]]$intensity, c(sumi(c(10, NA, 30)), sumi(c(NA, 25, 35))))
  expect_equal(res[[1]], res[[2]])
})

test_that(".process_peaks_data same-rt EIC cumsum with NAs", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200, 300), c(100, 200, 300), compress = FALSE),
    intensity = NumericList(c(10, NA, 30), c(NA, 25, 35), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(
    rtMin = c(1, 1), rtMax = c(2, 2),
    mzMin = c(100, 200), mzMax = c(200, 300)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  ## cumsum path strips NAs first, then sums within range
  ## spec1: non-NA mz=c(100,300) int=c(10,30) → mz 100-200: 10, mz 200-300: 30
  ## spec2: non-NA mz=c(200,300) int=c(25,35) → mz 100-200: 25, mz 200-300: 60
  expect_equal(res[[1]]$intensity, c(10, 25))
  expect_equal(res[[2]]$intensity, c(30, 60))
})

test_that(".process_peaks_data different-rt with NAs in intensities", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, NA), c(NA, 25), c(30, 40), compress = FALSE),
    rtime = c(1, 2, 3),
    msLevel = rep(1L, 3),
    dataOrigin = rep("A", 3)
  ))
  cd <- data.frame(
    rtMin = c(1, 2), rtMax = c(2, 3),
    mzMin = c(-Inf, -Inf), mzMax = c(Inf, Inf)
  )
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  ## TIC: sumi ignores NAs → spec1: 10, spec2: 25, spec3: 70
  expect_equal(res[[1]]$intensity, c(sumi(c(10, NA)), sumi(c(NA, 25))))
  expect_equal(res[[2]]$intensity, c(sumi(c(NA, 25)), sumi(c(30, 40))))
})

test_that(".process_peaks_data same-rt all-NA spectrum", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), compress = FALSE),
    intensity = NumericList(c(NA_real_, NA_real_), c(10, 20), compress = FALSE),
    rtime = c(1, 2),
    msLevel = rep(1L, 2),
    dataOrigin = rep("A", 2)
  ))
  cd <- data.frame(rtMin = 1, rtMax = 2, mzMin = 100, mzMax = 200)
  ## cumsum: spec1 all-NA → NA, spec2 → 30
  res <- .process_peaks_data(cd, sp, c("rtime", "intensity"), sumi, FALSE)
  expect_true(is.na(res[[1]]$intensity[1]))
  expect_equal(res[[1]]$intensity[2], 30)
  ## maxi: spec1 all-NA → NA, spec2 max = 20
  res_max <- .process_peaks_data(cd, sp, c("rtime", "intensity"), maxi, FALSE)
  expect_true(is.na(res_max[[1]]$intensity[1]))
  expect_equal(res_max[[1]]$intensity[2], 20)
})


test_that("intensity/rtime bypass peaksData when queue is empty", {
  ## Direct backend dispatch when no processing queue
  int_direct <- intensity(c_full)
  int_via_be <- peaksData(.backend(c_full), columns = "intensity", drop = TRUE)
  expect_identical(int_direct, int_via_be)

  rt_direct <- rtime(c_full)
  rt_via_be <- peaksData(.backend(c_full), columns = "rtime", drop = TRUE)
  expect_identical(rt_direct, rt_via_be)

  ## Same test with ChromBackendMzR
  int_mzr <- intensity(c_mzr)
  rt_mzr <- rtime(c_mzr)
  expect_identical(int_mzr,
    peaksData(.backend(c_mzr), columns = "intensity", drop = TRUE))
  expect_identical(rt_mzr,
    peaksData(.backend(c_mzr), columns = "rtime", drop = TRUE))

  ## With processing queue, should go through full peaksData dispatch
  c_queued <- filterPeaksData(c_full, variables = "rtime",
                              ranges = c(12.5, 45.5))
  int_queued <- intensity(c_queued)
  expect_false(identical(int_queued, int_direct))
})

test_that("lengths works for ChromBackendMemory", {
  l <- lengths(be)
  expect_equal(l, vapply(peaksData(be), nrow, integer(1L)))
  expect_type(l, "integer")
  ## Empty backend still has one chromatogram, with zero data points
  expect_equal(lengths(be_empty), 0L)
})

test_that("lengths works for Chromatograms", {
  ## No queue
  l_full <- lengths(c_full)
  expect_equal(l_full, vapply(peaksData(c_full), nrow, integer(1L)))

  l_mzr <- lengths(c_mzr)
  expect_equal(l_mzr, vapply(peaksData(c_mzr), nrow, integer(1L)))

  l_empty <- lengths(c_empty)
  expect_equal(l_empty, 0L)

  ## With queue
  c_queued <- filterPeaksData(c_full, variables = "rtime",
                              ranges = c(12.5, 14.0))
  l_queued <- lengths(c_queued)
  expect_equal(l_queued, vapply(peaksData(c_queued), nrow, integer(1L)))
})

test_that("peaksData fast [[ path for single column", {
  ## drop=TRUE, single column should use fast path
  int_fast <- peaksData(be, columns = "intensity", drop = TRUE)
  int_ref <- lapply(peaksData(be), `[[`, "intensity")
  expect_equal(int_fast, int_ref)

  rt_fast <- peaksData(be, columns = "rtime", drop = TRUE)
  rt_ref <- lapply(peaksData(be), `[[`, "rtime")
  expect_equal(rt_fast, rt_ref)
})

test_that("setBackend clears processing queue", {
  c_queued <- filterPeaksData(c_mzr, variables = "rtime",
                              ranges = c(12.5, 45.5))
  expect_length(.processingQueue(c_queued), 1L)
  pd_before <- peaksData(c_queued)

  c_mem <- setBackend(c_queued, ChromBackendMemory())
  expect_length(.processingQueue(c_mem), 0L)
  expect_identical(peaksData(c_mem), peaksData(.backend(c_mem)))
  expect_equal(peaksData(c_mem), pd_before)
})


test_that(".logging appends a timestamped entry", {
  res <- .logging(character(), "test message")
  expect_length(res, 1L)
  expect_match(res, "^test message \\[")

  res2 <- .logging(res, "second")
  expect_length(res2, 2L)
  expect_match(res2[2], "^second \\[")
})

test_that(".strip_na_peaks no NAs returns unchanged", {
  mz <- c(100, 200, 300)
  int <- c(10, 20, 30)
  res <- .strip_na_peaks(mz, int)
  expect_equal(res$mz, mz)
  expect_equal(res$int, int)
  expect_equal(res$n, 3L)
})

test_that(".strip_na_peaks partial NAs strips correctly", {
  mz <- c(100, 200, 300, 400)
  int <- c(10, NA, 30, NA)
  res <- .strip_na_peaks(mz, int)
  expect_equal(res$mz, c(100, 300))
  expect_equal(res$int, c(10, 30))
  expect_equal(res$n, 2L)
})

test_that(".strip_na_peaks all NAs returns NULL", {
  mz <- c(100, 200)
  int <- c(NA_real_, NA_real_)
  expect_null(.strip_na_peaks(mz, int))
})

test_that(".strip_na_peaks empty vectors returns unchanged", {
  res <- .strip_na_peaks(numeric(0), numeric(0))
  expect_equal(res$mz, numeric(0))
  expect_equal(res$int, numeric(0))
  expect_equal(res$n, 0L)
})

## ── .mz_boundaries ──────────────────────────────────────────────────────────

test_that(".mz_boundaries returns correct lo/hi indices", {
  mz <- c(100, 200, 300, 400, 500)

  ## Single range enclosing all
  b <- .mz_boundaries(100, 500, mz)
  expect_equal(b$lo, 1L)
  expect_equal(b$hi, 5L)

  ## Range in the middle
  b <- .mz_boundaries(200, 400, mz)
  expect_equal(b$lo, 2L)
  expect_equal(b$hi, 4L)

  ## Range excluding all (below)
  b <- .mz_boundaries(10, 50, mz)
  expect_true(b$lo > b$hi) # empty range

  ## Multiple ranges at once
  b <- .mz_boundaries(c(100, 300), c(200, 500), mz)
  expect_equal(b$lo, c(1L, 3L))
  expect_equal(b$hi, c(2L, 5L))
})

test_that(".rt_keep finite bounds returns correct indices", {
  rt <- c(1, 2, 3, 4, 5)
  res <- .rt_keep(rt, 2, 4, 5L)
  expect_equal(res, 2:4)
})

test_that(".rt_keep infinite bounds returns all", {
  rt <- c(1, 2, 3)
  res <- .rt_keep(rt, -Inf, Inf, 3L)
  expect_equal(res, 1:3)
})

test_that(".rt_keep narrow window returns subset", {
  rt <- c(1.0, 2.0, 3.0, 4.0)
  res <- .rt_keep(rt, 2.5, 3.5, 4L)
  expect_equal(res, 3L)
})

test_that(".rt_keep no match returns empty", {
  rt <- c(1, 2, 3)
  res <- .rt_keep(rt, 10, 20, 3L)
  expect_equal(res, integer(0))
})

## ── .make_chrom_entry ────────────────────────────────────────────────────────

test_that(".make_chrom_entry returns data.frame with both columns", {
  res <- .make_chrom_entry(
    rtime = c(1, 2, 3), intensity = c(10, 20, 30),
    do_rt = TRUE, do_int = TRUE, columns = c("rtime", "intensity"),
    drop = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(res$rtime, c(1, 2, 3))
  expect_equal(res$intensity, c(10, 20, 30))
  expect_equal(nrow(res), 3L)
})

test_that(".make_chrom_entry rtime-only", {
  res <- .make_chrom_entry(
    rtime = c(1, 2), intensity = numeric(0),
    do_rt = TRUE, do_int = FALSE, columns = "rtime", drop = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(names(res), "rtime")
})

test_that(".make_chrom_entry drop=TRUE single column returns vector", {
  res <- .make_chrom_entry(
    rtime = numeric(0), intensity = c(10, 20, 30),
    do_rt = FALSE, do_int = TRUE, columns = "intensity", drop = TRUE
  )
  expect_type(res, "double")
  expect_equal(res, c(10, 20, 30))
})

test_that(".make_chrom_entry drop=TRUE two columns returns data.frame", {
  res <- .make_chrom_entry(
    rtime = c(1, 2), intensity = c(10, 20),
    do_rt = TRUE, do_int = TRUE, columns = c("rtime", "intensity"),
    drop = TRUE
  )
  expect_s3_class(res, "data.frame")
})

## ── .spectra_format_chromData ────────────────────────────────────────────────

test_that(".spectra_format_chromData returns correct structure", {
  sp <- Spectra(DataFrame(
    mz = NumericList(c(100, 200), c(100, 200), compress = FALSE),
    intensity = NumericList(c(10, 20), c(15, 25), compress = FALSE),
    rtime = c(1.5, 3.0),
    msLevel = rep(1L, 2),
    dataOrigin = rep("fileA", 2),
    chromSpectraIndex = rep("grp1", 2)
  ))
  res <- .spectra_format_chromData(sp)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_equal(res$msLevel, 1L)
  expect_equal(res$rtMin, 1.5)
  expect_equal(res$rtMax, 3.0)
  expect_equal(res$mzMin, -Inf)
  expect_equal(res$mzMax, Inf)
  expect_equal(res$dataOrigin, "fileA")
  expect_equal(res$chromSpectraIndex, "grp1")
})

test_that(".peaksData errors for wrong object type", {
  expect_error(.peaksData("not_a_backend"),
               "must be of class 'Chromatograms' or 'ChromBackend'")
})

test_that(".chromData errors for wrong object type", {
  expect_error(.chromData(42),
               "must be of class 'Chromatograms' or 'ChromBackend'")
})

## ── .compare_chrom_pair ─────────────────────────────────────────────────────

test_that(".compare_chrom_pair returns 1 for identical chromatograms", {
    pd <- data.frame(rtime = c(1, 2, 3, 4), intensity = c(10, 50, 100, 50))
    expect_equal(.compare_chrom_pair(pd, pd), 1)
})

test_that(".compare_chrom_pair returns NA for non-overlapping RT", {
    pd_a <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    pd_b <- data.frame(rtime = c(10, 11, 12), intensity = c(40, 50, 60))
    expect_true(is.na(.compare_chrom_pair(pd_a, pd_b)))
})

test_that(".compare_chrom_pair returns NA when fewer than 2 points", {
    pd_short <- data.frame(rtime = 1, intensity = 50)
    pd_ok <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    expect_true(is.na(.compare_chrom_pair(pd_short, pd_ok)))
    expect_true(is.na(.compare_chrom_pair(pd_ok, pd_short)))
})

test_that(".compare_chrom_pair uses custom FUN when provided", {
    pd_a <- data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    pd_b <- data.frame(rtime = c(1, 2, 3), intensity = c(5, 15, 25))
    cosine <- function(x, y, ...) {
        sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    }
    res <- .compare_chrom_pair(pd_a, pd_b, FUN = cosine)
    expected <- cosine(pd_a$intensity, pd_b$intensity)
    expect_equal(res, expected, tolerance = 1e-10)
})

test_that(".compare_chrom_pair interpolates onto common RT grid", {
    ## Different RT grids, but overlapping range
    pd_a <- data.frame(rtime = c(1, 2, 3, 4), intensity = c(10, 20, 30, 40))
    pd_b <- data.frame(rtime = c(2, 3, 4, 5), intensity = c(20, 30, 40, 50))
    res <- .compare_chrom_pair(pd_a, pd_b)
    expect_true(is.numeric(res))
    expect_true(!is.na(res))
    ## Overlap is RT 2-4, perfectly correlated linear data → correlation = 1
    expect_equal(res, 1, tolerance = 1e-10)
})

test_that(".compare_chrom_pair respects method argument", {
    pd_a <- data.frame(rtime = c(1, 2, 3, 4, 5),
                       intensity = c(10, 20, 30, 20, 10))
    pd_b <- data.frame(rtime = c(1, 2, 3, 4, 5),
                       intensity = c(5, 15, 25, 15, 5))
    res_p <- .compare_chrom_pair(pd_a, pd_b, method = "pearson")
    res_s <- .compare_chrom_pair(pd_a, pd_b, method = "spearman")
    expect_true(is.numeric(res_p))
    expect_true(is.numeric(res_s))
    ## Both should be 1 for perfectly linearly related data
    expect_equal(res_p, 1, tolerance = 1e-10)
    expect_equal(res_s, 1, tolerance = 1e-10)
})

## ── .compare_chromatograms ──────────────────────────────────────────────────

test_that(".compare_chromatograms self-comparison returns symmetric matrix", {
    pd <- list(
        data.frame(rtime = c(1, 2, 3, 4), intensity = c(10, 50, 100, 50)),
        data.frame(rtime = c(1, 2, 3, 4), intensity = c(5, 25, 50, 25)),
        data.frame(rtime = c(10, 11, 12), intensity = c(40, 50, 60))
    )
    res <- .compare_chromatograms(pd)
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(3L, 3L))
    expect_equal(diag(res), c(1, 1, 1))
    expect_equal(res[1, 2], res[2, 1])
    expect_equal(res[1, 3], res[3, 1])
    expect_equal(res[2, 3], res[3, 2])
    ## Chromatograms 1 & 2 overlap and are correlated
    expect_true(!is.na(res[1, 2]))
    ## Chromatograms 1/2 vs 3 don't overlap
    expect_true(is.na(res[1, 3]))
    expect_true(is.na(res[2, 3]))
})

test_that(".compare_chromatograms cross-comparison returns n x m matrix", {
    pd_x <- list(
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30)),
        data.frame(rtime = c(1, 2, 3), intensity = c(5, 15, 25))
    )
    pd_y <- list(
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    )
    res <- .compare_chromatograms(pd_x, pd_y)
    expect_true(is.matrix(res))
    expect_equal(dim(res), c(2L, 1L))
    ## First pair is identical
    expect_equal(res[1, 1], 1, tolerance = 1e-10)
    ## Second pair is linearly related → correlation = 1
    expect_equal(res[2, 1], 1, tolerance = 1e-10)
})

test_that(".compare_chromatograms empty lists return 0-dim matrix", {
    expect_equal(dim(.compare_chromatograms(list())), c(0L, 0L))
    pd <- list(data.frame(rtime = 1:3, intensity = c(10, 20, 30)))
    res <- .compare_chromatograms(list(), pd)
    expect_equal(nrow(res), 0L)
    res2 <- .compare_chromatograms(pd, list())
    expect_equal(ncol(res2), 0L)
})

test_that(".compare_chromatograms applies labels for self-comparison", {
    pd <- list(
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30)),
        data.frame(rtime = c(1, 2, 3), intensity = c(5, 15, 25))
    )
    res <- .compare_chromatograms(pd, labels = c("a", "b"))
    expect_equal(rownames(res), c("a", "b"))
    expect_equal(colnames(res), c("a", "b"))
})

test_that(".compare_chromatograms with custom FUN", {
    cosine <- function(x, y, ...) {
        sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    }
    pd <- list(
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30)),
        data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30))
    )
    res <- .compare_chromatograms(pd, FUN = cosine)
    expect_equal(res[1, 2], 1, tolerance = 1e-10)
    expect_equal(res[2, 1], 1, tolerance = 1e-10)
})

test_that(".compare_chromatograms single element returns 1x1 with diag = 1", {
    pd <- list(data.frame(rtime = c(1, 2, 3), intensity = c(10, 20, 30)))
    res <- .compare_chromatograms(pd)
    expect_equal(dim(res), c(1L, 1L))
    expect_equal(res[1, 1], 1)
})

## ── .resolve_labels ─────────────────────────────────────────────────────────

test_that(".resolve_labels returns NULL when labels is NULL", {
    expect_null(.resolve_labels(c_full, NULL))
})

test_that(".resolve_labels returns values from a valid column", {
    res <- .resolve_labels(c_full, "mz")
    expect_equal(res, chromData(c_full)[["mz"]])
})

test_that(".resolve_labels errors on non-existent column", {
    expect_error(.resolve_labels(c_full, "nonexistent"), "not found")
})

test_that(".resolve_labels errors on duplicated column values", {
    c_dup <- c_full
    c_dup$msLevel <- rep(1L, length(c_dup))
    expect_error(.resolve_labels(c_dup, "msLevel"), "duplicated")
})

test_that(".resolve_labels errors on non-character or length != 1", {
    expect_error(.resolve_labels(c_full, 42), "single character string")
    expect_error(.resolve_labels(c_full, c("mz", "msLevel")),
                 "single character string")
})
