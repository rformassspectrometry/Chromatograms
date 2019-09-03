#' @include hidden_aliases.R
NULL

#' @description
#'
#' Check if chromData has all required columns.
#'
#' @noRd
#'
#' @param x chromData `DataFrame`
.valid_chrom_data_required_columns <- function(x, columns = c("dataStorage")) {
    if (nrow(x)) {
        missing_cn <- setdiff(columns, colnames(x))
        if (length(missing_cn))
            return(paste0("Required column(s): ",
                          paste(missing_cn, collapse = ", "),
                          " is/are missing"))
    }
    NULL
}

#' Function to check data types of selected columns in the provided `DataFrame`.
#'
#' @param x `DataFrame` to validate.
#'
#' @param datatypes named `character`, names being column names and elements
#'     expected data types.
#'
#' @author Johannes Rainer
#'
#' @noRd
.valid_column_datatype <- function(x, datatypes = .CHROMATOGRAMS_DATA_COLUMNS) {
    datatypes <- datatypes[names(datatypes) %in% colnames(x)]
    res <- mapply(FUN = .is_class, x[, names(datatypes), drop = FALSE],
                  datatypes)
    if (!all(res))
        paste0("The following columns have a wrong data type: ",
               paste(names(res[!res]), collapse = ", "),
               ". The expected data type(s) is/are: ",
               paste(datatypes[names(res)[!res]], collapse = ", "), ".")
    else NULL
}

#' @importFrom MsCoreUtils vapply1l
#'
#' @noRd
.valid_rtime_column <- function(x) {
    if (length(x$rtime)) {
        if (!all(vapply1l(x$rtime, is.numeric)))
            return("rtime column should contain a list of numeric")
        if (any(vapply1l(x$rtime, is.unsorted)))
            return("rtime values have to be sorted increasingly")
    }
    NULL
}

.valid_intensity_column <- function(x) {
    if (length(x$intensity))
        if (!all(vapply1l(x$intensity, is.numeric)))
            return("intensity column should contain a list of numeric")
    NULL
}

.valid_intensity_rtime_columns <- function(x) {
    ## Don't want to have that tested on all on-disk objects.
    if (length(x$intensity) && length(x$rtime))
        if (any(lengths(x$intensity) != lengths(x$rtime)))
            return(paste0("Length of rtime and intensity values differ for",
                          " some chromatograms"))
    NULL
}

.CHROMATOGRAMS_DATA_COLUMNS <- c(
    chromIndex = "integer",
    collisionEnergy = "numeric",
    dataOrigin = "character",
    dataStorage = "character",
    intensity = "NumericList",
    msLevel = "integer",
    mz = "numeric",
    mzMin = "numeric",
    mzMax = "numeric",
    precursorMz = "numeric",
    precursorMzMin = "numeric",
    precursorMzMax = "numeric",
    productMz = "numeric",
    productMzMin = "numeric",
    productMzMax = "numeric",
    rtime = "NumericList"
)

#' @rdname ChromBackend
#'
#' @export ChromBackendDataFrame
ChromBackendDataFrame <- function() {
    new("ChromBackendDataFrame")
}

#' Helper function to return a column from the (chrom data) `DataFrame`. If
#' the column `column` is an `Rle` `as.vector` is called on it. If column is
#' the name of a mandatory variable but it is not available it is created on
#' the fly.
#'
#' @note This function is equivalent to the `get_rle_column` function in
#'     the `Spectra` package (defined in *MsBackendDataFrame-functions.R*).
#'
#' @param x `DataFrame`
#'
#' @param column `character(1)` with the name of the column to return.
#'
#' @importMethodsFrom S4Vectors [[
#'
#' @author Johannes Rainer
#'
#' @noRd
.get_rle_column <- function(x, column) {
    if (any(colnames(x) == column)) {
        if (is(x[[column]], "Rle"))
            as.vector(x[[column]])
        else x[[column]]
    } else if (any(names(.CHROMATOGRAMS_DATA_COLUMNS) == column)) {
        nr_x <- nrow(x)
        if (nr_x)
            as(rep.int(NA, nr_x), .CHROMATOGRAMS_DATA_COLUMNS[column])
        else
            do.call(.CHROMATOGRAMS_DATA_COLUMNS[column], args = list())
    } else stop("column '", column, "' not available")
}

#' Helper function to combine backends that base on [ChromBackendDataFrame()].
#'
#' @param objects `list` of `ChromBackend` objects.
#'
#' @return [ChromBackend()] object with combined content.
#'
#' @author Johannes Rainer
#'
#' @importFrom MsCoreUtils rbindFill vapply1c
#'
#' @noRd
.combine_chrom_backend_data_frame <- function(objects) {
    if (length(objects) == 1)
        return(objects[[1]])
    if (!all(vapply1c(objects, class) == class(objects[[1]])))
        stop("Can only merge backends of the same type: ", class(objects[[1]]))
    res <- new(class(objects[[1]]))
    suppressWarnings(
        res@chromData <- asRleDataFrame(do.call(
            rbindFill, lapply(objects, function(z) z@chromData)),
            columns = c("dataStorage", "dataOrigin"))
    )
    if (any(colnames(res@chromData) == "rtime"))
        res@chromData$rtime[is.na(res@chromData$rtime)] <- list(numeric())
    if (any(colnames(res@chromData) == "intensity"))
        res@chromData$intensity[is.na(res@chromData$intensity)] <-
            list(numeric())
    res
}
