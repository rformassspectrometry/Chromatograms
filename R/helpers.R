## Here i want to group all the helper function used in the package...
# put explanation that this is for the reused fct and list methd in which it is used

#' @note
#' Used for:
#' - backendMerge for ChromBackendMemory, I actually do not know how to make it
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


#' Helper function that matches `x` against `mz` (using the `closest` function)
#' and returns the indices of `x` that match any of the values in `mz`. The
#' function takes care of sorting `x` and `mz` and deals also with missing
#' values.
#'
#' @return `integer` with the indices of values in `x` that are not `NA` and
#'     are matching any of the values in `mz` given `ppm` and `tolerance`.
#'
#' @noRd
#'
#' @importFrom MsCoreUtils common
#'
#' @note
#' Used in:
#' - filterMzValues
#'
#' @author Sebastian Gibb, Johannes Rainer
.values_match_mz <- function(x, mz, ppm = 20, tolerance = 0) {
    o <- order(x, na.last = NA)
    cmn <- common(x[o], sort(mz), tolerance = tolerance, ppm = ppm,
                  duplicates = "keep", .check = FALSE)
    sort(o[cmn])
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
#' - `ValidPeaksData()`
#' @noRd
.check_column_order_and_types <- function(df, expected_cols, expected_types) {
    # Check if the first two columns are 'rtime' and 'intensity'
    if (!identical(colnames(df)[1:2], expected_cols))
        return(paste0("Columns should be in the order 'rtime', 'intensity'."))
    # Check if the data types match the expected types
    invalid_cols <- sapply(expected_cols, function(col) {
        !is(df[[col]], expected_types[[col]])
    })
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
#' - `ValidPeaksData()`
#' @noRd
.check_rtime <- function(df) {
    if (nrow(df) == 0) return(NULL)
    if (any(is.na(df$rtime)))
        return("'rtime' column contains NA values.")

    if (!all(diff(df$rtime) > 0))
        return("'rtime' column is not strictly increasing.")

    return(NULL)
}

#' Function to validate each peaksData entry
#'
#' @note:
#' used in:
#' - `ValidPeaksData()`
#' @noRd
.validate_entry <- function(df, i, expected_cols, expected_types) {
    msgs <- NULL
    if (!is.data.frame(df))
        msgs <- c(msgs, paste0("Entry ", i, ": all 'peaksData' entries should ",
                               "be of class 'data.frame'"))
    else
        msgs <- c(msgs,
                  .check_column_order_and_types(df,
                                                expected_cols,
                                                expected_types),
                  .check_rtime(df))
    return(msgs)
}
