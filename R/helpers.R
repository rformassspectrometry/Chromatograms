#' Here are the helper functions used in the package.
#' Please add a description of the function and the methods in which it is used.

#' @note
#' Used for:
#' - `backendMerge()` for ChromBackendMemory, I actually do not know how to make it
#' so that it applies to other future backend. I guess this is something to
#' think about.
#'
#' @author Johannes Rainer
#' @importFrom MsCoreUtils vapply1c rbindFill
#' @noRd
.df_combine <- function(objects) {
    if (length(objects) == 1)
        return(objects[[1]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    res <- objects[[1]]
    pv <- names(res@peaksData[[1]])
    for (i in 2:length(objects)) {
        res@chromData <- rbindFill(res@chromData, objects[[i]]@chromData)
        pv2 <- peaksVariables(objects[[i]])
        if (length(pv) == length(pv2) && all(pv == pv2)) {
            res@peaksData <- c(res@peaksData , objects[[i]]@peaksData)
        } else
            stop("Provided objects have different sets of peak variables. ",
                 "Combining such objects is currently not supported.")
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
    if (!identical(colnames(df)[1:2], expected_cols))
        return(paste0("Columns should be in the order 'rtime', 'intensity'."))
    invalid_cols <- vapply(expected_cols, function(col) {
        !is(df[[col]], expected_types[[col]]) }, logical(1))
    if (any(invalid_cols)) {
        invalid_col_names <- expected_cols[invalid_cols]
        return(paste0("The peaksData variable(s) ", paste(invalid_col_names,
                                                          collapse = ", "),
                      " have the wrong data type."))
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
    if (nrow(df) == 0) return(NULL)
    if (any(is.na(df$rtime)))
        return("'rtime' column contains NA values.")

    if (!all(diff(df$rtime) > 0))
        return("'rtime' column is not strictly increasing.")

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
        for (i in seq_along(queue))
            object <- do.call(queue[[i]]@FUN, c(object, queue[[i]]@ARGS))
        return(object)
    }
    if (!is(f, "factor")) stop("f must be a factor")
    if (length(f) != length(object))
        stop("length 'f' has to be equal to the length of 'object' (",
             length(object), ")")
    processed_data <- bplapply(split(object, f), function(x) {
        for (i in seq_along(queue))
            x <- do.call(queue[[i]]@FUN, c(x, queue[[i]]@ARGS))
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
    if (!is.data.frame(df))
        msgs <- c(msgs, paste0("Entry ", i, ": all 'peaksData' entries should ",
                               "be of class 'data.frame'"))
    else
        msgs <- c(msgs, .check_column_order_and_types(df, expected_cols,
                                                      expected_types),
                  .check_rtime(df))
    return(msgs)
}

#' Function to validate the processingQueue slot of a Chromatograms object
#'
#' Used in:
#' - `validObject(Chromatograms())`
#' @importFrom MsCoreUtils vapply1l
#' @noRd
.valid_processing_queue <- function(x) {
    if (length(x) && !all(vapply1l(x, inherits, "ProcessingStep")))
        stop("'processingQueue' should only contain ProcessingStep objects.")
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
    if (length(ranges) != 2 * nc)
        stop("Length of 'ranges' needs to be twice the length of the ",
             "parameter 'query'")

    # Compute within_ranges for each column of the query
    within_ranges <- vapply(seq_len(nc), function(i) {
        pairs <- c(ranges[2 * i - 1], ranges[2 * i])
        between(query[[i]], pairs)
    }, logical(nrow(query)))

    if (match == "all") {
        if (nr == 1) return(as.integer(all(within_ranges)))
        return(which(rowSums(within_ranges) == nc))
    }
    if (nr == 1) return(as.integer(any(within_ranges)))
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
        if (!requireNamespace("mzR", quietly = TRUE))
            stop("The use of 'ChromBackendMzR' requires package 'mzR'. ",
                 "Install it using 'BiocManager::install(\"mzR\")'")
}

#' Used In:
#' - `ChromBackendSpectra()`
#' @noRd
.check_Spectra_package <- function() {
    if (!requireNamespace("Spectra", quietly = TRUE))
        stop("The use of 'ChromBackendSpectra' requires package 'Spectra'. ",
             "Install it using 'BiocManager::install(\"Spectra\")'")
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
    colnames(tmp)[colnames(tmp) == "chromatogramIndex"] <- "chromIndex"
    colnames(tmp)[colnames(tmp) == "precursorCollisionEnergy"] <- "collisionEnergy"
    colnames(tmp)[colnames(tmp) == "productIsolationWindowTargetMZ"] <- "productMz"
    colnames(tmp)[colnames(tmp) == "precursorIsolationWindowTargetMZ"] <- "precursorMz"
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
.plot_single_chromatogram <- function(x, xlab = "rtime (s)", ylab = "intensity",
                                      type = "l", xlim = numeric(),
                                      ylim = numeric(),
                                      main = paste("m/z", round(mz(x), 1)),
                                      col = "#00000080",add = FALSE,
                                      axes = TRUE, frame.plot = axes,
                                      orientation = 1, ...) {
    v <- peaksData(x)[[1L]]
    rts <- v$rtime
    ints <- orientation * v[, "intensity"]
    if (!length(xlim))
        suppressWarnings(xlim <- range(rts, na.rm = TRUE))
    if (!length(ylim))
        suppressWarnings(
            ylim <- range(orientation * c(0, max(abs(ints), na.rm = TRUE))))
    if (any(is.infinite(xlim)))
        xlim <- c(0, 0)
    if (any(is.infinite(ylim)))
        ylim <- c(0, 0)
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
        if (frame.plot)
            box(...)
        title(main = main, xlab = xlab, ylab = ylab, ...)
    }
    plot.xy(xy.coords(rts, ints), type = type, col = col, ...)
}

#' Used In:
#' - `peaksData` for `ChromBackendSpectra` class.
#' @noRd
.process_peaks_data <- function(cd, s, columns, fun, drop) {
    s <- Spectra::filterRanges(s, spectraVariables = rep("rtime", nrow(cd)),
                      ranges = as.vector(rbind(cd$rtmin, cd$rtmax)),
                      match = "any")
    pd <- Spectra::peaksData(s, columns = c("mz", "intensity"))
    do_rt <- "rtime" %in% columns
    do_int <- "intensity" %in% columns
    rt <- rtime(s)
    lapply(seq_len(nrow(cd)), function(i) {
        keep <- between(rt, c(cd$rtmin[i], cd$rtmax[i]))
        df <- as.data.frame(matrix(ncol = 0, nrow = sum(keep)))
        if (do_rt)
            df$rtime <- rt[keep]
        if (do_int)
            df$intensity <- vapply(pd[keep], function(z) {
                fun(z[between(z[, "mz"], c(cd$mzMin[i], cd$mzMax[i])),
                      "intensity"])
            }, NA_real_, USE.NAMES = FALSE)
        df[, columns, drop = drop]
    })
}

#' Used in:
#' - `backendInitialize()` for `ChrombackendSpectra`
#' @noRd
.spectra_format_chromData <- function(sps) {
    data.frame(
        msLevel = unique(sps$msLevel),
        rtmin = min(sps$rtime, na.rm = TRUE),
        rtmax = max(sps$rtime, na.rm = TRUE),
        mzMin = -Inf,
        mzMax = Inf,
        mz = Inf,
        polarity = sps$polarity[1],
        scanWindowLowerLimit = sps$scanWindowLowerLimit[1],
        scanWindowUpperLimit = sps$scanWindowUpperLimit[1],
        dataOrigin = unique(sps$dataOrigin),
        chromSpectraIndex = unique(sps$chromSpectraIndex)
    )
}

#' Used in:
#' - `factorize()` for `ChrombackendSpectra`
#' @noRd
.ensure_rt_mz_columns <- function(chrom_data, spectra, spectra_f) {
    ## Ensure mzmin and mzmax are either both present or both missing
    if (!all(c("mzmin", "mzmax") %in% colnames(chrom_data))) {
        if ("mzmin" %in% colnames(chrom_data) || "mzmax" %in% colnames(chrom_data)) {
            stop("Both 'mzmin' and 'mzmax' must be present if one is provided.")
        } else {
            chrom_data$mzmin <- -Inf
            chrom_data$mzmax <- Inf
        }
    }

    ## Ensure rtmin and rtmax are either both present or computed
    if (!all(c("rtmin", "rtmax") %in% colnames(chrom_data))) {
        if ("rtmin" %in% colnames(chrom_data) || "rtmax" %in% colnames(chrom_data)) {
            stop("Both 'rtmin' and 'rtmax' must be present if one is provided.")
        } else {
            rt_range <- lapply(split(spectra$rtime, spectra_f), function(rt) {
                list(rtmin = min(rt, na.rm = TRUE), rtmax = max(rt, na.rm = TRUE))
            })
            rt_values <- do.call(rbind, rt_range)
            chrom_data$rtmin <- rt_values[, "rtmin"]
            chrom_data$rtmax <- rt_values[, "rtmax"]
        }
    }
    chrom_data
}
