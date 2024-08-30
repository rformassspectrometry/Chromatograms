#' @rdname ChromBackend
#'
#' @export
coreChromVariables <- function() .CORE_CHROM_VARIABLES

#' *core* chromatogram variables with expected data type: `integer`, `numeric`,
#' and `character`. Must be a single value.
#'
#' @noRd
.CORE_CHROM_VARIABLES <- c(
    chromIndex = "integer",
    collisionEnergy = "numeric",
    dataOrigin = "character",
    dataStorage = "character",
    msLevel = "integer",
    mz = "numeric",
    mzMin = "numeric",
    mzMax = "numeric",
    precursorMz = "numeric",
    precursorMzMin = "numeric",
    precursorMzMax = "numeric",
    productMz = "numeric",
    productMzMin = "numeric",
    productMzMax = "numeric"
    )

#' @title Fill data frame with columns for missing core chrom variables
#'
#' @description
#'
#' `fillCoreChromVariables()` fills a provided `data.frame`
#' with columns for eventually missing *core* chromatogram variables.
#' The missing core variables are added as new columns with missing values
#' (`NA`) of the correct data type.
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
    cv <- .CORE_CHROM_VARIABLES
    miss <- cv[setdiff(names(cv), colnames(x))]
    if (!length(miss))
        return(x)
    cbind(x, lapply(miss, function(z, n) rep(as(NA, z), n), nr))
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
#' @importFrom methods is
#'
#' @export
validChromData <- function(x = data.frame(), error = TRUE) {
    .valid_chrom_backend_data_storage(x$dataStorage)
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

#' *core* peaks variables with expected data type:`numeric`, and
#' must be plural.
#'
#' @noRd
.CORE_PEAKS_VARIABLES <- c(
    intensity = "numeric",
    rtime = "numeric"
)

#' @rdname ChromBackend
#'
#' @export
corePeaksVariables <- function() .CORE_PEAKS_VARIABLES

#' `validPeaksData()` checks that the names of the input peaksData list,
#' representing *core* peaks variables are of the correct data type.
#'
#' @param x `list` representing the peaks data of a `Chromatograms`
#'
#' @importFrom methods is
#'
#' @export

#' @rdname hidden_aliases
validPeaksData <- function(x = list()) {
    if (!is.list(x)) {
        stop("'peaksData' must be a 'list'")
    }
    lapply(x, function(df) {
        if (!is.data.frame(df))
            stop("All the peaksData entries should be of class 'data.frame'")
        pv <- intersect(colnames(df), names(.CORE_PEAKS_VARIABLES))
        if (!all(vapply(pv, function(col) is(df[[col]]
                                             , .CORE_PEAKS_VARIABLES[[col]]),
                        logical(1))))
            stop("The peaksData variable ", col, " has the wrong data type.")
    })
}

#' @noRd
.valid_chrom_backend_data_storage <- function(x) {
    if (anyNA(x))
        return("'NA' values in dataStorage are not allowed.")
    NULL
}  ## what's this ?
