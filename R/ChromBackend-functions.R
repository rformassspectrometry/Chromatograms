.valid_chrom_backend_data_storage <- function(x) {
    if (anyNA(x))
        return("'NA' values in dataStorage are not allowed.")
    NULL
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
#' @author Sebastian Gibb, Johannes Rainer
.values_match_mz <- function(x, mz, ppm = 20, tolerance = 0) {
    o <- order(x, na.last = NA)
    cmn <- common(x[o], sort(mz), tolerance = tolerance, ppm = ppm,
                  duplicates = "keep", .check = FALSE)
    sort(o[cmn])
}

#' @rdname ChromBackend
#'
#' @export
coreChromVariables <- function() .CORE_CHROM_VARIABLES

#' *core* chromatogram variables with expected data type:
#'
#' @noRd
.CORE_CHROM_VARIABLES <- c(
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

#' @title Fill data frame with columns for missing core chrom variables
#'
#' @description
#'
#' `fillCoreChromVariables()` fills a provided `data.frame` (or `DataFrame`)
#' with columns for eventually missing *core* chromatogram variables **except**
#' peaks variables (i.e. `"intensity"` and `"rtime"`). The missing core
#' variables are added as new columns with missing values (`NA`) of the
#' correct data type.
#' Use [coreChromVariables()] to list the set of core variables and their data
#' types.
#'
#' @param x `data.frame` with potentially present core chrom variable columns
#'
#' @return input data frame `x` with missing core variables added (with the
#'     correct data type).
#'
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#'
#' ## Define a data frame
#' a <- data.frame(msLevel = c(1L, 1L, 2L), other_column = "b")
#'
#' ## Add missing core chromatogram variables to this data frame
#' fillCoreChromVariables(a)
#'
#' ## The data.frame thus contains columns for all core chromatogram
#' ## variables in the respective expected data type (but filled with
#' ## missing values).
fillCoreChromVariables <- function(x = data.frame()) {
    nr <- nrow(x)
    cv <- .CORE_CHROM_VARIABLES[!names(.CORE_CHROM_VARIABLES) %in%
                                c("intensity", "rtime")]
    miss <- setdiff(names(cv), colnames(x))
    cbind(x, lapply(cv, function(z, n) rep(as(NA, z), n), nr))
}

#' @title Check core chromatogram variables for correct data types
#'
#' @description
#'
#' `validChromData()` checks that columns, representing *core* chromatogram
#' variables are of the correct data type.
#'
#' @param x `data.frame` representing metadata of a `Chromatograms`
#'
#' @param error `logical(1)` whether an error should be thrown (the default)
#'     if one or more columns don't have the correct data type.
#'
#' @return
#'
#' If core variables have all the correct data type: an empty character.
#' If one or more core variables (columns) have the wrong data type the
#' function either throws an error (with `error = TRUE`) or returns a
#' `character` specifying which variables/columns don't have the correct
#' type (for `error = FALSE`).
#'
#' @export
validChromData <- function(x = data.frame(), error = TRUE) {
    cn <- intersect(colnames(x), names(.CORE_CHROM_VARIABLES))
    msg <- unlist(lapply(cn, function(z) {
        if (!is(x[, z], .CORE_CHROM_VARIABLES[z]))
            paste0("Column \"", z, "\" has the wrong data type. ")
        else character()
    }), use.names = FALSE)
    if (length(msg) && error)
        stop(msg)
    else msg
}

#' Create a list of empty peak matrices
#'
#' @param x `integer` with the number of (empty) matrices to create.
#'
#' @param columns `character` with the column names for each peak matrix
#'
#' @return `list` of length `x` with 0-row peak matrices
#'
#' @noRd
.empty_peaks_data <- function(x, columns = c("mz", "intensity")) {
    emat <- matrix(numeric(), ncol = length(columns), nrow = 0,
                   dimnames = list(character(), columns))
    replicate(x, emat, simplify = FALSE)
}
