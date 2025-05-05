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

#' @title Fill data.frame with columns for missing core chromatogram variables.
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
#'
#' @rdname hidden_aliases
fillCoreChromVariables <- function(x = data.frame()) {
    nr <- nrow(x)
    cv <- .CORE_CHROM_VARIABLES
    miss <- cv[setdiff(names(cv), colnames(x))]
    if (!length(miss)) {
        return(x)
    }
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
#' If the core variables have all the correct data type: an empty character.
#' If one or more core variables (columns) have the wrong data type the
#' function either throws an error (with `error = TRUE`) or returns a
#' `character` specifying which variables/columns don't have the correct
#' type (for `error = FALSE`).
#'
#' @importFrom methods is
#'
#' @export
#' @rdname hidden_aliases
validChromData <- function(x = data.frame(), error = TRUE) {
    cn <- intersect(colnames(x), names(.CORE_CHROM_VARIABLES))
    msg <- unlist(lapply(cn, function(z) {
        if (!is(x[, z], .CORE_CHROM_VARIABLES[z])) {
            paste0("Column \"", z, "\" has the wrong data type. ")
        } else {
            NULL
        }
    }), use.names = FALSE)
    if (length(msg) && error) {
        stop(msg)
    } else {
        msg
    }
}

#' *core* peaks variables with expected data type:`numeric`, and
#' must be plural.
#'
#' @noRd
.CORE_PEAKS_VARIABLES <- c(
    rtime = "numeric",
    intensity = "numeric"
)

#' An empty peaks data.frame
#' @noRd
.EMPTY_PEAKS_DATA <- as.data.frame(lapply(
    .CORE_PEAKS_VARIABLES,
    function(x) vector(x, 0)
))

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
validPeaksData <- function(x = list(), error = TRUE) {
    if (!is.list(x)) stop("'peaksData' must be a 'list'")
    if (!length(x)) {
        return(NULL)
    }
    first_cols <- colnames(x[[1]])
    expected_cols <- names(.CORE_PEAKS_VARIABLES)
    expected_types <- .CORE_PEAKS_VARIABLES
    msgs <- unlist(lapply(seq_along(x), function(i) {
        df <- x[[i]]
        if (!identical(colnames(df), first_cols)) {
            return(paste(
                "All data.frames must have the same columns in the",
                " same order. Issue found in entry", i
            ))
        }
        .validate_entry(x[[i]], i, expected_cols, expected_types)
    }))
    if (length(msgs) && error) {
        stop(msgs)
    } else {
        msgs
    }
}
