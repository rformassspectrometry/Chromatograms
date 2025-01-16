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

#' Helper function that checks if the `dataStorage` of a `Chromatogram` object
#' contains any `NA` values.
#'
#' @note
#' Used in:
#' - `validChromData()`
#' @noRd
.valid_chrom_backend_data_storage <- function(x) {
    if (anyNA(x))
        return("'NA' values in dataStorage are not allowed.")
    NULL
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
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
.run_process_queue <- function(object,
                               f = processingChunkFactor(object),
                               BPPARAM = SerialParam()) {
    queue <- object@processingQueue
    if (!length(queue)) return(object@backend)
    BPPARAM <- backendBpparam(object@backend, BPPARAM)
    if (!length(f) || length(levels(f)) == 1) {
        for (i in seq_along(queue))
            object@backend <- do.call(queue[[i]]@FUN, c(object@backend, queue[[i]]@ARGS))
        return(object@backend )
    }
    if (!is(f, "factor")) stop("f must be a factor")
    if (length(f) != length(object))
        stop("length 'f' has to be equal to the length of 'object' (",
             length(object), ")")
    processed_data <- bplapply(split(object@backend, f), function(x) {
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
        return("'processingQueue' should only contain ProcessingStep objects.")
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
