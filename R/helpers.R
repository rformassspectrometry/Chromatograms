#' Here are the helper functions used in the package.
#' Please add a description of the function and the methods in which it is
#' used.

#' @note
#' Used for:
#' - `backendMerge()`
#'
#' @author Johannes Rainer
#' @importFrom MsCoreUtils vapply1c rbindFill
#' @noRd
.df_combine <- function(objects) {
    if (length(objects) == 1) {
        return(objects[[1]])
    }
    if (!all(vapply1c(objects, class) == class(objects[[1]]))) {
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    }
    res <- objects[[1]]
    pv <- names(.peaksData(res)[[1]])
    for (i in 2:length(objects)) {
        res@chromData <- rbindFill(.chromData(res), .chromData(objects[[i]]))
        pv2 <- peaksVariables(objects[[i]])
        if (length(pv) == length(pv2) && all(pv == pv2)) {
            res@peaksData <- c(.peaksData(res), .peaksData(objects[[i]]))
        } else {
            stop(
                "Provided objects have different sets of peak variables. ",
                "Combining such objects is currently not supported."
            )
        }
    }
    res
}

#' Fast rbind for lists of data.frames using `data.table::rbindlist()`.
#' Always returns a plain `data.frame`.
#' @importFrom data.table rbindlist
#' @noRd
.fast_rbind <- function(lst) {
  if (length(lst) == 0L) return(data.frame())
  as.data.frame(rbindlist(lst, use.names = TRUE, fill = TRUE))
}

#' Helper function to check the order and data types of columns
#'
#' @note:
#' used in:
#' - `validPeaksData()`
#' @noRd
.check_column_order_and_types <- function(df, expected_cols, expected_types) {
    if (!identical(colnames(df)[seq_len(2)], expected_cols)) {
        return(paste0("Columns should be in the order 'rtime', 'intensity'."))
    }
    invalid_cols <- vapply(expected_cols, function(col) {
        !is(df[[col]], expected_types[[col]])
    }, logical(1))
    if (any(invalid_cols)) {
        invalid_col_names <- expected_cols[invalid_cols]
        return(paste0(
            "The peaksData variable(s) ", paste(invalid_col_names,
                collapse = ", "
            ),
            " have the wrong data type."
        ))
    }
    return(NULL)
}

#' Helper function to check the properties of the 'rtime' column.
#'
#' @note:
#' used in:
#' - `validPeaksData()`
#' @noRd
.check_rtime <- function(df) {
    if (nrow(df) == 0) {
        return(NULL)
    }
    if (any(is.na(df$rtime))) {
        return("'rtime' column contains NA values.")
    }

    if (!all(diff(df$rtime) > 0)) {
        return("'rtime' column is not strictly increasing.") ## does it need to strictly increase ?
    }

    return(NULL)
}

#' Function to apply the processing queue to the backend, return a peaksData.
#' Used in:
#' - `peaksData(Chromatograms())`
#' - `applyProcessing()`
#'
#' It takes a backend and a preprocessingQueue and applies it. It returns
#' then the backend. This function might need to be  refined later in case the
#' backend is `readOnly == TRUE`.
#'
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
.run_process_queue <- function(object, queue, f = factor(),
    BPPARAM = SerialParam()) {
    BPPARAM <- backendBpparam(object, BPPARAM)
    if (!length(f) || length(levels(f)) == 1) {
        for (i in seq_along(queue)) {
            object <- do.call(queue[[i]]@FUN, c(object, queue[[i]]@ARGS))
        }
        return(object)
    }
    if (!is(f, "factor")) stop("f must be a factor")
    if (length(f) != length(object)) {
        stop(
            "length 'f' has to be equal to the length of 'object' (",
            length(object), ")"
        )
    }
    processed_data <- bplapply(split(object, f), function(x) {
        for (i in seq_along(queue)) {
            x <- do.call(queue[[i]]@FUN, c(x, queue[[i]]@ARGS))
        }
        x
    }, BPPARAM = BPPARAM)
    backendMerge(processed_data)
}

#' Function to validate each peaksData entry
#'
#' @note:
#' used in:
#' - `validPeaksData()`
#' @noRd
.validate_entry <- function(df, i, expected_cols, expected_types) {
    msgs <- NULL
    if (!is.data.frame(df)) {
        msgs <- c(msgs, paste0(
            "Entry ", i, ": all 'peaksData' ",
            "entries should ",
            "be of class 'data.frame'"
        ))
    } else {
        msgs <- c(
            msgs, .check_column_order_and_types(
                df, expected_cols,
                expected_types
            ),
            .check_rtime(df)
        )
    }
    return(msgs)
}

#' Function to validate the processingQueue slot of a `Chromatograms` object
#'
#' Used in:
#' - `validObject(Chromatograms())`
#' @importFrom MsCoreUtils vapply1l
#' @noRd
.valid_processing_queue <- function(x) {
    if (length(x) && !all(vapply1l(x, inherits, "ProcessingStep"))) {
        stop("'processingQueue' should only contain ProcessingStep objects.")
    }
    NULL
}

#' function to loop through  query column and check if within corresponding
#' ranges. Return an index of the corresponding matches.
#' Used in:
#' - `filterPeaksData()`: looped through the list of data.frame
#' - `filterChromData()`
#' @importFrom MsCoreUtils between
#' @noRd
.filter_ranges <- function(query, ranges, match) {
    nc <- ncol(query)
    nr <- nrow(query)
    if (length(ranges) != 2 * nc) {
        stop(
            "Length of 'ranges' needs to be twice the length of the ",
            "parameter 'query'"
        )
    }

    # Compute within_ranges for each column of the query
    within_ranges <- vapply(seq_len(nc), function(i) {
        pairs <- c(ranges[2 * i - 1], ranges[2 * i])
        between(query[[i]], pairs)
    }, logical(nrow(query)))

    if (match == "all") {
        if (nr == 1) {
            return(as.integer(all(within_ranges)))
        }
        return(which(rowSums(within_ranges) == nc))
    }
    if (nr == 1) {
        return(as.integer(any(within_ranges)))
    }
    return(which(rowSums(within_ranges) > 0))
}

#' Used in:
#' - `filterPeaksData()`
#' @noRd
.logging <- function(x, ...) {
    c(x, paste0(..., " [", date(), "]"))
}

#' Used In:
#' - `ChromBackendMzR()`
#' @noRd
.check_mzR_package <- function() {
    if (!requireNamespace("mzR", quietly = TRUE)) {
        stop(
            "The use of 'ChromBackendMzR' requires package 'mzR'. ",
            "Install it using 'BiocManager::install(\"mzR\")'"
        )
    }
}

#' Function to create chromData form mzml file
#' Used In:
#' - `backendInitialize()` for `ChromBackendMzR` class
#' Helper function to format chromatographic data from mzR files.
#' @noRd
.mzR_format_chromData <- function(file) {
    .check_mzR_package()
    msd <- mzR::openMSfile(file)
    on.exit(mzR::close(msd))
    tmp <- mzR::chromatogramHeader(msd)
    colnames(tmp)[colnames(tmp) ==
        "chromatogramIndex"] <- "chromIndex"
    colnames(tmp)[colnames(tmp) ==
        "precursorCollisionEnergy"] <- "collisionEnergy"
    colnames(tmp)[colnames(tmp) ==
        "productIsolationWindowTargetMZ"] <- "productMz"
    colnames(tmp)[colnames(tmp) ==
        "precursorIsolationWindowTargetMZ"] <- "precursorMz"
    tmp$dataOrigin <- file
    tmp
}

#' Used In:
#' - `peaksData()` for `ChromBackendMzR` class
#' @noRd
.get_chrom_data <- function(fl, idx) {
    .check_mzR_package()
    msd <- mzR::openMSfile(fl)
    on.exit(mzR::close(msd))
    mzR::chromatogram(msd, idx, drop = FALSE)
}

#' Helper function to plot a single chromatogram.
#' @note:
#' Used in:
#' - `plotChromatograms()`
#' - `plotChromatogramsOverlay()`
#'
#' @importFrom graphics plot.new plot.window plot.xy axis box title par
#' @importFrom grDevices dev.hold dev.flush xy.coords n2mfrow
#' @noRd
.plot_single_chromatogram <- function(x, xlab = "rtime (s)",
    ylab = "intensity",
    type = "l", xlim = numeric(),
    ylim = numeric(),
    main = paste("m/z", round(mz(x), 1)),
    col = "#00000080", add = FALSE,
    axes = TRUE, frame.plot = axes,
    orientation = 1, ...) {
    v <- peaksData(x)[[1L]]
    rts <- v$rtime
    raw_ints <- v[, "intensity"]
    ints <- orientation * raw_ints
    if (!length(xlim)) {
        xlim <- range(rts, na.rm = TRUE)
    }
    if (!length(ylim)) {
        ylim <- range(orientation * c(0, max(abs(ints), na.rm = TRUE)))
    }
    if (any(is.infinite(xlim))) {
        xlim <- c(0, 0)
    }
    if (any(is.infinite(ylim))) {
        ylim <- c(0, 0)
    }
    if (!add) {
        dev.hold()
        on.exit(dev.flush())
        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
    }
    if (!add) {
        if (axes) {
            axis(side = 1, ...)
            axis(side = 2, ...)
        }
        if (frame.plot) {
            box(...)
        }
        title(main = main, xlab = xlab, ylab = ylab, ...)
    }
    plot.xy(xy.coords(rts, ints), type = type, col = col, ...)
}
#' Strip NA intensities from paired mz/intensity vectors.
#'
#' Returns a list with cleaned `mz`, `int`, and `n`, or `NULL` if no
#' non-NA values remain.
#'
#' Used in:
#' - `.build_intensity_matrix()`
#' - `.compute_chrom_intensities()`
#' @noRd
.strip_na_peaks <- function(mz_j, int_j) {
  if (!anyNA(int_j))
    return(list(mz = mz_j, int = int_j, n = length(mz_j)))
  keep <- !is.na(int_j)
  int_j <- int_j[keep]
  mz_j <- mz_j[keep]
  n <- length(int_j)
  if (n == 0L) NULL else list(mz = mz_j, int = int_j, n = n)
}

#' Compute m/z range boundary indices via binary search.
#'
#' Wraps `findInterval()` to return the lo (first included) and hi (last
#' included) positions for each m/z range in `mzMins`/`mzMaxs`.
#'
#' Used in:
#' - `.build_intensity_matrix()`
#' - `.compute_chrom_intensities()`
#' @noRd
.mz_boundaries <- function(mzMins, mzMaxs, mz_j) {
  list(
    lo = findInterval(mzMins, mz_j, left.open = TRUE,
                      checkSorted = FALSE, checkNA = FALSE) + 1L,
    hi = findInterval(mzMaxs, mz_j,
                      checkSorted = FALSE, checkNA = FALSE)
  )
}

#' Find spectra indices within a retention time window.
#'
#' Used in:
#' - `.process_peaks_data()`
#' @noRd
.rt_keep <- function(valid_rt, rt_min, rt_max, n_valid) {
  if (is.finite(rt_min) || is.finite(rt_max))
    which(valid_rt >= rt_min & valid_rt <= rt_max)
  else
    seq_len(n_valid)
}

#' Build a single chromatogram result entry.
#'
#' Constructs a data.frame (or vector when `drop = TRUE`) directly without
#' calling `as.data.frame()`, which avoids its per-call overhead (~10x
#' faster for repeated construction).
#'
#' Used in:
#' - `.process_peaks_data()`
#' @noRd
.make_chrom_entry <- function(rtime, intensity, do_rt, do_int, columns, drop) {
  df <- list()
  if (do_rt) df$rtime <- rtime
  if (do_int) df$intensity <- intensity
  class(df) <- "data.frame"
  attr(df, "row.names") <- .set_row_names(
    if (do_rt) length(rtime) else length(intensity)
  )
  if (drop && length(columns) == 1L) df[[columns]] else df
}

#' Pre-extract and pre-filter spectra data for chromatogram computation.
#'
#' Extracts mz, intensity, and rtime from a Spectra object as plain R lists,
#' removes duplicated retention times, and pre-filters spectra to the global
#' retention time range. Also classifies the request (shared rt, TIC/BPC,
#' cumsum-eligible) to guide downstream dispatch.
#'
#' Used in:
#' - `.process_peaks_data()`
#'
#' @param cd `data.frame` with columns rtMin, rtMax, mzMin, mzMax.
#' @param s `Spectra` object.
#' @param fun Summary function (e.g., `sumi`, `maxi`).
#' @return A named `list` with prepared vectors and classification flags.
#' @importFrom Spectra mz intensity rtime
#' @importFrom MsCoreUtils sumi maxi
#' @noRd
.prepare_spectra_input <- function(cd, s, fun) {
  ## as.list() converts SimpleNumericList to plain list (~12x faster for [[)
  mz_list <- as.list(mz(s))
  int_list <- as.list(intensity(s))
  rt <- rtime(s)
  not_dup <- !duplicated(rt)
  global_rt_min <- min(cd$rtMin)
  global_rt_max <- max(cd$rtMax)
  global_keep <- not_dup & rt >= global_rt_min & rt <= global_rt_max
  valid_idx <- which(global_keep)
  rt_mins <- cd$rtMin
  rt_maxs <- cd$rtMax
  mzMins <- cd$mzMin
  mzMaxs <- cd$mzMax
  list(
    valid_rt   = rt[valid_idx],
    valid_mz   = mz_list[valid_idx],
    valid_int  = int_list[valid_idx],
    n_valid    = length(valid_idx),
    rt_mins    = rt_mins,
    rt_maxs    = rt_maxs,
    mzMins     = mzMins,
    mzMaxs     = mzMaxs,
    n_cd       = nrow(cd),
    same_rt    = all(rt_mins == rt_mins[1L]) && all(rt_maxs == rt_maxs[1L]),
    all_tic    = all(is.infinite(mzMins) & mzMins < 0) &&
                 all(is.infinite(mzMaxs) & mzMaxs > 0),
    use_cumsum = identical(fun, sumi)
  )
}

#' Compute intensities for one spectrum across one or more m/z ranges.
#'
#' Core per-spectrum logic shared by `.build_intensity_matrix()` (matrix
#' path) and `.compute_chrom_intensities()` (vector path).  Dispatches
#' through three branches: TIC/BPC (no m/z filtering), cumsum (for
#' `sumi`), and generic (e.g. `maxi`).
#'
#' Used in:
#' - `.build_intensity_matrix()`
#' - `.compute_chrom_intensities()`
#'
#' @param mz_j,int_j Numeric vectors of m/z and intensity for one spectrum.
#' @param mzMins,mzMaxs Numeric vectors of m/z bounds (length = n ranges).
#' @param fun Summary function.
#' @param all_tic Logical: all mz ranges infinite (TIC/BPC case)?
#' @param use_cumsum Logical: use cumsum optimisation (for `sumi`)?
#' @return Numeric vector of length `length(mzMins)`.
#' @noRd
.spectrum_intensities <- function(mz_j, int_j, mzMins, mzMaxs,
                                  fun, all_tic, use_cumsum) {
  n_cd <- length(mzMins)
  if (all_tic) return(rep(fun(int_j), n_cd))
  if (use_cumsum) {
    pk <- .strip_na_peaks(mz_j, int_j)
    if (is.null(pk)) return(rep(NA_real_, n_cd))
    cs <- c(0, cumsum(pk$int))
    b <- .mz_boundaries(mzMins, mzMaxs, pk$mz)
    res <- rep(NA_real_, n_cd)
    w <- which(b$lo <= b$hi & b$lo >= 1L & b$hi >= 1L & b$lo <= pk$n)
    if (length(w)) res[w] <- cs[b$hi[w] + 1L] - cs[b$lo[w]]
    return(res)
  }
  n_mz <- length(mz_j)
  if (n_mz == 0L) return(rep(NA_real_, n_cd))
  b <- .mz_boundaries(mzMins, mzMaxs, mz_j)
  vapply(seq_len(n_cd), function(i) {
    if (b$lo[i] <= b$hi[i] && b$lo[i] >= 1L && b$hi[i] <= n_mz)
      fun(int_j[b$lo[i]:b$hi[i]])
    else NA_real_
  }, numeric(1))
}

#' Build intensity matrix for the shared-retention-time path.
#'
#' When all chromatograms share the same retention time window, computes
#' intensities using a spectrum-major loop via `.spectrum_intensities()`.
#'
#' Used in:
#' - `.process_peaks_data()`
#'
#' @param valid_mz,valid_int Lists of m/z and intensity vectors.
#' @param kept Integer indices into `valid_*` for spectra in the rt window.
#' @param mzMins,mzMaxs Numeric vectors of m/z bounds per chromatogram.
#' @param n_cd Number of chromatograms.
#' @param fun Summary function.
#' @param all_tic Logical: all mz ranges infinite (TIC/BPC case)?
#' @param use_cumsum Logical: use cumsum optimisation (for `sumi`)?
#' @return A `matrix` of dimensions `length(kept)` x `n_cd`.
#' @noRd
.build_intensity_matrix <- function(valid_mz, valid_int, kept, mzMins,
                                    mzMaxs, n_cd, fun, all_tic, use_cumsum) {
  n_kept <- length(kept)
  int_mat <- matrix(NA_real_, nrow = n_kept, ncol = n_cd)
  for (j in seq_len(n_kept))
    int_mat[j, ] <- .spectrum_intensities(
      valid_mz[[kept[j]]], valid_int[[kept[j]]],
      mzMins, mzMaxs, fun, all_tic, use_cumsum)
  int_mat
}

#' Compute intensity vector for a single chromatogram (different-rt path).
#'
#' Thin wrapper around `.spectrum_intensities()` that returns a vector
#' (one value per spectrum) for a single m/z range.
#'
#' Used in:
#' - `.process_peaks_data()`
#'
#' @param valid_mz,valid_int Lists of m/z and intensity vectors.
#' @param kept Integer indices into `valid_*` for spectra in this chrom's
#'   rt window.
#' @param mz_lo,mz_hi Scalar m/z bounds for this chromatogram.
#' @param fun Summary function.
#' @param use_cumsum Logical: use cumsum optimisation?
#' @return Numeric vector of length `length(kept)`.
#' @noRd
.compute_chrom_intensities <- function(valid_mz, valid_int, kept, mz_lo,
                                       mz_hi, fun, use_cumsum) {
  all_tic <- is.infinite(mz_lo) && mz_lo < 0 && is.infinite(mz_hi)
  vapply(kept, function(j)
    .spectrum_intensities(valid_mz[[j]], valid_int[[j]],
                          mz_lo, mz_hi, fun, all_tic, use_cumsum),
    numeric(1), USE.NAMES = FALSE)
}

#' Compute peaks data from Spectra for ChromBackendSpectra.
#'
#' Main orchestrator that converts raw spectra into chromatographic peaks data
#' by aggregating intensities within defined m/z and retention time ranges.
#' Dispatches to `.build_intensity_matrix()` (shared-rt path) or
#' `.compute_chrom_intensities()` (per-chromatogram path).
#'
#' Used In:
#' - `peaksData` for `ChromBackendSpectra` class.
#' @importFrom Spectra peaksData
#' @importFrom MsCoreUtils maxi
#' @noRd
.process_peaks_data <- function(cd, s, columns, fun, drop) {
  do_rt <- "rtime" %in% columns
  do_int <- "intensity" %in% columns
  prep <- .prepare_spectra_input(cd, s, fun)
  n_cd <- prep$n_cd
  valid_rt <- prep$valid_rt
  valid_mz <- prep$valid_mz
  valid_int <- prep$valid_int
  n_valid <- prep$n_valid
  if (prep$same_rt) {
    ## Optimized path: all chromatograms share the same rt range
    kept <- .rt_keep(valid_rt, prep$rt_mins[1L], prep$rt_maxs[1L], n_valid)
    n_kept <- length(kept)
    rtime_out <- valid_rt[kept]
    int_mat <- if (do_int && n_kept > 0L)
      .build_intensity_matrix(valid_mz, valid_int, kept, prep$mzMins,
                              prep$mzMaxs, n_cd, fun, prep$all_tic,
                              prep$use_cumsum)
    res <- lapply(seq_len(n_cd), function(i)
      .make_chrom_entry(
        rtime_out,
        if (do_int && n_kept > 0L) int_mat[, i] else numeric(0L),
        do_rt, do_int, columns, drop
      ))
  } else {
    ## Fallback: per-chromatogram rt ranges differ
    res <- vector("list", n_cd)
    for (i in seq_len(n_cd)) {
      kept <- .rt_keep(valid_rt, prep$rt_mins[i], prep$rt_maxs[i], n_valid)
      n_kept <- length(kept)
      rtime_out <- valid_rt[kept]
      intensities <- if (do_int && n_kept > 0L)
        .compute_chrom_intensities(valid_mz, valid_int, kept,
                                   prep$mzMins[i], prep$mzMaxs[i],
                                   fun, prep$use_cumsum)
      else numeric(n_kept)
      res[[i]] <- .make_chrom_entry(
        rtime_out, intensities, do_rt, do_int, columns, drop
      )
    }
  }
  res
}

#' Used in:
#' - `backendInitialize()` for `ChrombackendSpectra`
#' @noRd
.spectra_format_chromData <- function(sps) {
    res <- data.frame(
        msLevel = unique(sps$msLevel),
        rtMin = min(sps$rtime, na.rm = TRUE),
        rtMax = max(sps$rtime, na.rm = TRUE),
        mzMin = -Inf,
        mzMax = Inf,
        mz = Inf,
        dataOrigin = unique(sps$dataOrigin),
        chromSpectraIndex = unique(sps$chromSpectraIndex)
    )
    ## Add optional columns if present
    if ("polarity" %in% spectraVariables(sps)) {
        res$polarity <- sps$polarity[1]
    }
    if ("scanWindowLowerLimit" %in% spectraVariables(sps)) {
        res$scanWindowLowerLimit <- sps$scanWindowLowerLimit[1]
    }
    if ("scanWindowUpperLimit" %in% spectraVariables(sps)) {
        res$scanWindowUpperLimit <- sps$scanWindowUpperLimit[1]
    }
    res
}

#' Used in:
#' - `factorize()` for `ChrombackendSpectra`
#' @noRd
.ensure_rt_mz_columns <- function(chrom_data, spectra, spectra_f) {
    if (!all(c("mzMin", "mzMax") %in% colnames(chrom_data))) {
        if ("mzMin" %in% colnames(chrom_data) ||
            "mzMax" %in% colnames(chrom_data)) {
            stop("Both 'mzMin' and 'mzMax' must be present if one",
                 " is provided.")
        } else {
            chrom_data$mzMin <- -Inf
            chrom_data$mzMax <- Inf
        }
    }
    if (!all(c("rtMin", "rtMax") %in% colnames(chrom_data))) {
        if ("rtMin" %in% colnames(chrom_data) || "rtMax" %in%
            colnames(chrom_data)) {
            stop("Both 'rtMin' and 'rtMax' must be present if one",
                 " is provided.")
        } else {
            levs <- levels(spectra_f)
            if (is.null(levs)) {
                levs <- unique(as.character(spectra_f))
            }
            rt_mat <- vapply(levs, function(lvl) {
                range(spectra$rtime[spectra_f == lvl], na.rm = TRUE)
            }, numeric(2))
            chrom_idx <- as.character(chrom_data$chromSpectraIndex)
            chrom_data$rtMin <- rt_mat[1, chrom_idx]
            chrom_data$rtMax <- rt_mat[2, chrom_idx]
        }
    }
    chrom_data
}

#' Used in:
#' - `chromExtract()`.
#' @noRd
.validate_chromExtract_input <- function(object,
                                         peak.table,
                                         by,
                                         required_cols = c("rtMin", "rtMax",
                                                           by)) {
    cd <- .chromData(object)
    if (!all(required_cols %in% names(peak.table))) {
        stop("`peak.table` must contain columns: ", paste(required_cols,
                                                          collapse = ", "), ".")
    }

    if (anyNA(peak.table$rtMin) || anyNA(peak.table$rtMax)) {
        stop("Columns 'rtMin' and 'rtMax' in `peak.table` cannot ",
             "contain NA values.")
    }
    if (!all(by %in% names(cd))) {
        stop("All 'by' columns must be present in chromData(object).")
    }
    if (nrow(cd) != nrow(unique(cd[, by, drop = FALSE]))) {
        stop("Combinations of 'by' columns must uniquely identify rows ",
             "in chromData.")
    }
}


#' Used in:
#' - `chromExract()`
#' @noRd
.match_chromdata_peaktable <- function(object, peak.table, by) {
    cd <- .chromData(object)
    chrom_keys <- interaction(cd[, by, drop = FALSE], drop = TRUE)
    peak_keys  <- interaction(peak.table[, by, drop = FALSE], drop = TRUE)

    # ensure all peak.table keys exist in chromData
    missing_keys <- setdiff(levels(peak_keys), levels(chrom_keys))
    if (length(missing_keys)) {
        stop("Some combinations in `peak.table` do not exist in chromData: ",
             paste(missing_keys, collapse = ", "))
    }

    ## Subset chromdata and only keep the row of interest.
    keep_idx <- chrom_keys %in% peak_keys
    object <- object[keep_idx]
    chrom_keys <- droplevels(chrom_keys[keep_idx])

    # align factor levels (so splitting matches between cd and peak.table)
    shared_levels <- intersect(levels(peak_keys), levels(chrom_keys))
    chrom_keys <- factor(as.character(chrom_keys), levels = shared_levels)
    peak_keys  <- factor(as.character(peak_keys),  levels = shared_levels)

    list(object = object, chrom_keys = chrom_keys, peak_keys = peak_keys)
}

#' Used in:
#' - `chromExtract()`
#' @noRd
.check_overl_columns <- function(object, peak.table, required_cols) {
    overl_cols <- names(peak.table) %in% chromVariables(object)
    extra_cols <- setdiff(names(peak.table)[overl_cols], required_cols)
    if (length(extra_cols)) {
        warning( "The following columns in `peak.table` already exist in ",
                 "`chromData` and will be replaced in the output: ",
                 paste(extra_cols, collapse = ", ")
        )
    }
    overl_cols
}


#' Used in:
#' - `imputePeaksData()`
#' @importFrom stats approx filter loess spline dnorm sd predict
#' @noRd
.impute <- function(x, method, window = 2, span = 0.25, sd = 1,
                    extrapolate = FALSE) {
    if (all(is.na(x))) return(x)

    na_idx <- which(is.na(x))
    if (length(na_idx) == 0) return(x)
    not_na_idx <- which(!is.na(x))
    x_out <- seq_along(x)
    # rule = 2 extrapolates, rule = 1 returns NA outside range
    approx_rule <- if (extrapolate) 2L else 1L
    x[na_idx] <- switch(method,
        linear = approx(not_na_idx, x[not_na_idx], xout = na_idx,
                        rule = approx_rule)$y,
        spline = {
            vals <- spline(not_na_idx, x[not_na_idx],
                           xout = na_idx, method = "natural")$y
            if (!extrapolate) {
                # Set values outside data range to NA
                outside <- na_idx < min(not_na_idx) | na_idx > max(not_na_idx)
                vals[outside] <- NA
            }
            vals
        },
        gaussian = {
            # Create symmetric Gaussian kernel
            kernel_range <- -window:window
            w <- dnorm(kernel_range, mean = 0, sd = sd)
            w <- w / sum(w)

            # Fill missing with linear approx to allow smoothing
            x_filled <- x
            x_filled[is.na(x_filled)] <- approx(not_na_idx, x[not_na_idx],
                                                xout = which(is.na(x_filled)),
                                                rule = approx_rule)$y
            smoothed <- filter(x_filled, filter = w, sides = 2,
                               circular = FALSE)
            smoothed[na_idx]
        },
        loess = {
            fit <- loess(x[not_na_idx] ~ not_na_idx, span = span)
            vals <- predict(fit, newdata = na_idx)
            if (!extrapolate) {
                outside <- na_idx < min(not_na_idx) | na_idx > max(not_na_idx)
                vals[outside] <- NA
            }
            vals
        }
    )
    # Fallback for any remaining NAs (only if extrapolate = TRUE)
    na_remaining <- is.na(x)
    if (extrapolate && any(na_remaining)) {
        warning("Method chosen could not fill all NAs. ",
                "Falling back to linear interpolation ",
                "for these positions.")
        x[na_remaining] <- approx(not_na_idx, x[not_na_idx],
                                         xout = which(na_remaining),
                                         rule = 2)$y
    }
    x
}


## Used in:
## - BackendInitialize, ChromBackendSpectra method
#' @noRd
.map_spectra_vars <- function(object, spectraVariables) {
    # Check variable validity for mapping from the Spectra object stored in the internal @spectra slot
    spectra <- .spectra(object)
    cd <- .chromData(object)
    if (!all(spectraVariables %in% spectraVariables(spectra)))
        stop("All 'spectraVariables' must exist in the Spectra object.")
    if (any(spectraVariables %in% colnames(cd))) {
        existing <- intersect(spectraVariables, colnames(cd))
        non_replaceable <- vapply(existing, function(v) !all(is.na(cd[[v]])),
                                  logical(1))
        if (any(non_replaceable)) {
            stop("None of the 'spectraVariables' must already exist in the ",
                 "chromData data.frame.")
        }
    }
    idx <- spectra$chromSpectraIndex
    spd <- spectraData(spectra, columns = spectraVariables)

    ## Aggregate and simplify singletons
    aggregated <- as.data.frame(
        lapply(spectraVariables, function(var) {
            res <- tapply(spd[[var]], idx, unique, simplify = FALSE)
            if (all(lengths(res) == 1L)) {
                res <- unlist(res, use.names = TRUE)
            }
            res
        }),
        stringsAsFactors = FALSE
    )
    names(aggregated) <- spectraVariables

    ## match order and combine
    aggregated <- aggregated[as.character(cd$chromSpectraIndex), , drop = FALSE]
    cd <- cbind(cd, aggregated)
    rownames(cd) <- NULL
    object@chromData <- cd
    object
}


#' Compute pairwise similarity between two chromatograms.
#'
#' Interpolates both chromatograms onto the union of their retention time
#' points within the overlapping range, then computes the correlation (or a
#' custom similarity via `FUN`).
#'
#' Used in:
#' - `compareChromatograms()`
#'
#' @param x,y `data.frame` with columns `rtime` and `intensity`.
#' @param method `character(1)` correlation method passed to `cor()`.
#' @param MAPFUN function to align two chromatograms' retention times.
#' @param FUN function to compute similarity from aligned intensities.
#' @param ... additional arguments passed to both `MAPFUN` and `FUN`.
#' @return `numeric(1)` similarity value.
#' @noRd
.compare_chrom_pair <- function(x, y, MAPFUN = matchRtime, FUN = cor, ...) {
    aligned <- MAPFUN(x, y, ...)
    if (length(aligned$x) < 2L) return(NA_real_)
    FUN(aligned$x, aligned$y, ...)
}

#' Compute a pairwise similarity matrix between two lists of peaks
#' data.frames.
#'
#' Compares each chromatogram in `pd_x` with each chromatogram in `pd_y`.
#' When `pd_y` is not provided (default), `pd_x` is compared against itself
#' and the result is a symmetric n x n matrix.
#'
#' Used in:
#' - `compareChromatograms()`
#'
#' @param pd_x `list` of `data.frame` each with columns `rtime` and
#'        `intensity`.
#' @param pd_y `list` of `data.frame` each with columns `rtime` and
#'        `intensity`. Defaults to `pd_x` (self-comparison).
#' @param MAPFUN function to align retention times.
#' @param FUN function to compute similarity.
#' @param labels optional `character` vector of row/column names. Only
#'        meaningful for self-comparison (when `pd_y` is not supplied).
#' @param ... passed to `.compare_chrom_pair()`.
#' @return A numeric `matrix` with `length(pd_x)` rows and `length(pd_y)`
#'         columns.
#' @noRd
.compare_chromatograms <- function(pd_x, pd_y = pd_x,
                                   MAPFUN = matchRtime, FUN = cor,
                                   labels = NULL, ...) {
    nx <- length(pd_x)
    ny <- length(pd_y)
    if (nx == 0L || ny == 0L)
        return(matrix(numeric(0), nx, ny))
    self <- identical(pd_x, pd_y)
    mat <- matrix(NA_real_, nx, ny)
    if (self) {
        diag(mat) <- 1
        for (i in seq_len(nx - 1L)) {
            for (j in (i + 1L):ny) {
                val <- .compare_chrom_pair(
                    pd_x[[i]], pd_y[[j]], MAPFUN = MAPFUN, FUN = FUN, ...)
                mat[i, j] <- val
                mat[j, i] <- val
            }
        }
    } else {
        for (i in seq_len(nx))
            for (j in seq_len(ny))
                mat[i, j] <- .compare_chrom_pair(
                    pd_x[[i]], pd_y[[j]], MAPFUN = MAPFUN, FUN = FUN, ...)
    }
    if (!is.null(labels)) {
        rownames(mat) <- labels
        colnames(mat) <- labels
    }
    mat
}

#' Resolve and validate the labels vector from a chromData column.
#'
#' Used in:
#' - `compareChromatograms()`
#'
#' @param object A `Chromatograms` object.
#' @param labels `character(1)` column name in `chromData()`, or `NULL`.
#' @return A character vector of labels, or `NULL`.
#' @noRd
.resolve_labels <- function(object, labels) {
    if (is.null(labels)) return(NULL)
    if (!is.character(labels) || length(labels) != 1L)
        stop("'labels' must be a single character string")
    labs <- chromData(object)[[labels]]
    if (is.null(labs))
        stop("Column '", labels, "' not found in chromData")
    if (anyDuplicated(labs))
        stop("Column '", labels, "' contains duplicated values")
    labs
}


#' Below are internal accessor functions, used ubiquitously in the package.
#' These directly access the internal slots (e.g., `@chromData`, `@peaksData`,
#' `@spectra`, etc.).
#' These are NOT to be used by general users.
#' @noRd
.backend <- function(object) {
    object@backend
}
.peaksData <- function(object) {
    if (is(object, "Chromatograms")) {
        return(object@backend@peaksData)
    }
    if (is(object, "ChromBackend")) {
        return(object@peaksData)
    }
    stop("'object' must be of class 'Chromatograms' or 'ChromBackend'.")
}
.chromData <- function(object) {
    if (is(object, "Chromatograms")) {
        return(object@backend@chromData)
    }
    if (is(object, "ChromBackend")) {
        return(object@chromData)
    }
    stop("'object' must be of class 'Chromatograms' or 'ChromBackend'.")
}
.inMemory <- function(object) {
    object@inMemory
}
.processing <- function(object) {
    object@processing
}
.processingQueue <- function(object) {
    object@processingQueue
}
.spectra <- function(object) {
    object@spectra
}

#' Compute peak boundaries for a single chromatogram.
#'
#' Finds the retention time boundaries of the tallest peak in a chromatogram
#' using `MsCoreUtils::valleys()` to locate the flanking valleys. When the
#' valley intensities exceed a baseline-relative threshold, falls back to a
#' threshold-based boundary search.
#'
#' Used in:
#' - `peakBoundary()`
#'
#' @param rtime numeric vector of retention times.
#' @param intensity numeric vector of intensities.
#' @param threshold numeric(1), fraction of the peak height above baseline
#'        used as a fallback cut-off (default 0.1).
#' @param baselineThreshold numeric(1), fraction of the peak height above
#'        baseline used to validate valley positions (default 0.1).
#' @param baselineQuantile numeric(1), quantile of the intensity distribution
#'        used as the baseline estimate (default 0.1).
#' @return A numeric vector of length 2 (left boundary, right boundary
#'   retention times), or `c(NA_real_, NA_real_)` when no valid peak
#'   is found.
#' @importFrom MsCoreUtils valleys
#' @importFrom stats quantile
#' @noRd
.peak_boundary_one <- function(rtime, intensity, threshold = 0.1,
                               baselineThreshold = 0.1,
                               baselineQuantile = 0.1) {
    n <- length(intensity)
    na_res <- c(NA_real_, NA_real_)
    if (n < 3L || all(is.na(intensity)))
        return(na_res)
    max_int <- max(intensity, na.rm = TRUE)
    if (max_int == 0)
        return(na_res)
    max_idx <- which.max(intensity)
    baseline_int <- quantile(intensity, probs = baselineQuantile,
                             na.rm = TRUE)
    peak_height <- max_int - baseline_int
    baseline_thresh <- baseline_int + peak_height * baselineThreshold
    v <- valleys(intensity, max_idx)
    left_idx  <- if ("left" %in% colnames(v)) v[1L, "left"] else 1L
    right_idx <- if ("right" %in% colnames(v)) v[1L, "right"] else n
    left_ok <- !is.na(intensity[left_idx]) &&
        intensity[left_idx] <= baseline_thresh &&
        !(left_idx > 1L && is.na(intensity[left_idx - 1L]))
    right_ok <- !is.na(intensity[right_idx]) &&
        intensity[right_idx] <= baseline_thresh &&
        !(right_idx < n && is.na(intensity[right_idx + 1L]))
    if (!left_ok || !right_ok) {
        thresh_val <- baseline_int + peak_height * threshold
        left_cand  <- which(intensity[seq_len(max_idx)] <= thresh_val)
        right_cand <- which(intensity[max_idx:n] <= thresh_val)
        left_idx  <- if (length(left_cand)) max(left_cand) else 1L
        right_idx <- if (length(right_cand))
            max_idx + min(right_cand) - 1L else n
    }
    c(unname(rtime[left_idx]), unname(rtime[right_idx]))
}
