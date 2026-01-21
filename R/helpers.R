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

#' Function to validate the processingQueue slot of a Chromatograms object
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
#' Used In:
#' - `peaksData` for `ChromBackendSpectra` class.
#' @importFrom Spectra peaksData filterRanges
#' @noRd
.process_peaks_data <- function(cd, s, columns, fun, drop) {
    ## Handle single spectrum case: filterRanges fails with length(s) == 1
    if (length(s) > 1) {
        s <- filterRanges(s,
            spectraVariables = rep("rtime", nrow(cd)),
            ranges = as.vector(rbind(cd$rtMin, cd$rtMax)),
            match = "any"
        )
    } else {
        ## For single spectrum, manually filter by rtime range
        if (length(s) == 1) {
            rt_in_range <- s$rtime >= min(cd$rtMin) & s$rtime <= max(cd$rtMax)
            if (!rt_in_range) {
                s <- s[integer(0)]  ## Return empty Spectra
            }
        }
    }
    pd <- peaksData(s, columns = c("mz", "intensity"))
    do_rt <- "rtime" %in% columns
    do_int <- "intensity" %in% columns
    rt <- rtime(s)
    lapply(seq_len(nrow(cd)), function(i) {
        ## only keep the first rt if there is duplication
        keep <- between(rt, c(cd$rtMin[i], cd$rtMax[i])) & !duplicated(rt)
        df <- as.data.frame(matrix(ncol = 0, nrow = sum(keep)))
        if (do_rt) {
            df$rtime <- rt[keep]
        }
        if (do_int) {
            df$intensity <- vapply(pd[keep], function(z) {
                fun(z[
                    between(z[, "mz"], c(cd$mzMin[i], cd$mzMax[i])),
                    "intensity"
                ])
            }, NA_real_, USE.NAMES = FALSE)
        }
        df[, columns, drop = drop]
    })
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
.impute <- function(x, method,
                    window = 2, span = 0.25, sd = 1) {
    if (all(is.na(x))) return(x)

    na_idx <- which(is.na(x))
    if (length(na_idx) == 0) return(x)

    not_na_idx <- which(!is.na(x))
    x_out <- seq_along(x)

    x[na_idx] <- switch(method,
        linear = approx(not_na_idx, x[not_na_idx],
                        xout = na_idx, rule = 2)$y,

        spline = spline(not_na_idx, x[not_na_idx],
                        xout = na_idx, method = "natural")$y,
        gaussian = {
            # Create symmetric Gaussian kernel
            kernel_range <- -window:window
            w <- dnorm(kernel_range, mean = 0, sd = sd)
            w <- w / sum(w)

            # Fill missing with linear approx to allow smoothing
            x_filled <- x
            x_filled[is.na(x_filled)] <- approx(not_na_idx, x[not_na_idx],
                                                xout = which(is.na(x_filled)),
                                                rule = 2)$y
            smoothed <- filter(x_filled, filter = w, sides = 2,
                               circular = FALSE)
            smoothed[na_idx]
        },
        loess = {
            fit <- loess(x[not_na_idx] ~ not_na_idx, span = span)
            predict(fit, newdata = na_idx)
        }
    )
    # Fallback for any remaining NAs
    na_remaining <- is.na(x)
    if (any(na_remaining)) {
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
## - BackendInitialize, chrombackendSPectra method
#' @noRd
.map_spectra_vars <- function(object, spectraVariables) {
    ## check variable validity
    spectra <- .spectra(object)
    cd <- .chromData(object)
    if (!all(spectraVariables %in% spectraVariables(spectra)))
        stop("All 'spectraVariables' must exist in 'spectra'.")
    if (any(spectraVariables %in% colnames(cd))) {
        existing <- intersect(spectraVariables, colnames(cd))
        non_replaceable <- vapply(existing, function(v) !all(is.na(cd[[v]])), logical(1))
        if (any(non_replaceable)) {
            stop("None of the 'spectraVariables' must already exist in 'chromData'.")
        }
    }
    idx <- spectra$chromSpectraIndex
    spd <- spectraData(spectra, columns = spectraVariables)

    ## Aggregate and simplify singletons
    aggregated <- as.data.frame(
        lapply(spectraVariables, function(var) {
            res <- tapply(spd[[var]], idx, unique, simplify = FALSE)
            ## If each element is length 1, unlist to atomic vector
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


## Below are internal accessors functions, these are used ubiquitously in the
## package. They directly access the slots. these are NOT to be used by general
## users.
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
