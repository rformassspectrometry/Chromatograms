#' @include helpers.R
#' @include hidden_aliases.R
#' @include ChromBackend.R
NULL

#' @title Chromatographic Data Backend for Spectra Objects
#'
#' @name ChromBackendSpectra
#'
#' @description
#' The `ChromBackendSpectra` class extends `ChromBackendMemory`, inheriting
#' all its slots and methods while providing additional functionality for
#' summarizing chromatographic data from [Spectra::Spectra()] objects.
#'
#' It can be initialized with a `Spectra` object, which is stored in the
#' `spectra` slot of the backend. Users can also provide a `data.frame`
#' containing chromatographic metadata, stored in `chromData`. This metadata
#' filters the `Spectra` object and generates `peaksData`. If `chromData` is
#' not provided, a default `data.frame` is created from the `Spectra` data.
#'
#' The *dataOrigin* core chromatogram variable should reflect the *dataOrigin*
#' of the `Spectra` object. The `factorize.by` parameter defines the variables
#' for grouping `Spectra` data into chromatographic data. The default is
#' `c("msLevel", "dataOrigin")`, which will defines separate chromatograms for
#' each combination of `msLevel` and `dataOrigin`. These variables must be in
#' both `Spectra` and `chromData` (if provided).
#'
#' The `summarize.method` parameter defines how spectral data intensity is
#' summarized:
#' - **"sum"**: Sums intensity to create a Total Ion Chromatogram (TIC).
#' - **"max"**: Takes max intensity for a Base Peak Chromatogram (BPC).
#'
#' If `chromData` or its factorization columns are modified, the `factorize()`
#' method must be called to update `chromSpectraIndex`.
#'
#' @details
#' No `peaksData` is stored until the user calls a function that generates it
#' (e.g., `rtime()`, `peaksData()`, `intensity()`). The `peaksData` slot
#' replacement is unsupported—modifications are temporary to optimize memory.
#' The `inMemory` slot indicates this with `TRUE`.
#'
#' `ChromBackendSpectra` should reuse `ChromBackendMemory` methods whenever
#' possible to keep implementations simple.
#'
#' @param chromData A `data.frame` with chromatographic data for use in
#'        `backendInitialize()`. If missing, a default is generated. Columns
#'        like `rtmin`, `rtmax`, `mzmin`, and `mzmax` must be provided and not
#'        contain `NA` values. Use `-Inf/Inf` for unspecified values. The
#'        `"dataOrigin"` column must match the `Spectra` object's `"dataOrigin"`.
#'
#' @param factorize.by A `character` vector of variables for grouping `Spectra`
#'        data into chromatographic data. Default: `c("msLevel", "dataOrigin")`.
#'        If `chromData` is provided, it must contain these columns.
#'
#' @param object A `ChromBackendSpectra` object.
#'
#' @param spectra A `Spectra` object.
#'
#' @param summarize.method A `character` string specifying intensity summary:
#'        `"sum"` (default) or `"max"`.
#'
#' @param ... Additional parameters.
#'
#' @author Philippine Louail, Johannes Rainer.
#'
#' @exportClass ChromBackendSpectra
NULL


#' @noRd
ChromBackendSpectra <- setClass(
    "ChromBackendSpectra",
    contains = "ChromBackendMemory",
    slots = c(
        inMemory = "logical",
        spectra = "Spectra",
        summaryFun = "function"
    ),
    prototype = prototype(
        chromData = fillCoreChromVariables(data.frame()),
        peaksData = list(.EMPTY_PEAKS_DATA),
        readonly = TRUE,
        spectra = Spectra::Spectra(),
        version = "0.1",
        inMemory = FALSE,
        summaryFun = sumi
    )
)

#' @rdname ChromBackendSpectra
#' @importFrom methods new
#' @export ChromBackendSpectra
ChromBackendSpectra <- function() {
    .check_Spectra_package()
    new("ChromBackendSpectra")
}

#' @rdname ChromBackendSpectra
#' @importFrom methods callNextMethod
#' @importFrom MsCoreUtils rbindFill sumi maxi
setMethod("backendInitialize", "ChromBackendSpectra",
          function(object, spectra = Spectra::Spectra(),
                   factorize.by = c("msLevel" , "dataOrigin"),
                   summarize.method = c("sum", "max"),
                   chromData = data.frame(),
                   ...) {
              summarize.method <- match.arg(summarize.method)
              object@summaryFun <- if (summarize.method == "sum") sumi else maxi
              if (!is(spectra, "Spectra"))
                  stop("'spectra' must be a 'Spectra' object.")
              if (!length(spectra)) return(object)
              if (!all(factorize.by %in% Spectra::spectraVariables(spectra)))
                  stop("All 'factorize.by' variables must exist in 'spectra'.")
              if (!is.data.frame(chromData))
                  stop("'chromData' must be a 'data.frame'.")

              if (nrow(chromData) > 0) {
                  if (!all(factorize.by %in% colnames(chromData)))
                      stop("All 'factorize.by' variables must exist in ",
                           "'chromData'.")
                  object <- callNextMethod(object, chromData = chromData)
                  object <- factorize(object, factorize.by)
                  if (all(is.na(object$chromIndex)))
                    object$chromIndex <- seq_len(nrow(chromData))
                  object@spectra <- spectra
                  return(object)
                  #I have an issue here if the user provide a chromData but
                  # no rtmin/rtmax, these are not corevariables so they are not
                  # automatically added. and for mzmin mzmax evenif the columns
                  # are automaticall created,
                  # the code later WILL FAIL if the input is NA.
                  # Not sure how to handle that..
              }

              f <- factor(
                  do.call(
                      paste,
                      c(as.list(Spectra::spectraData(spectra)[, factorize.by]),
                        sep = "_")))
              spectra$chromSpectraIndex <- f
              full_sp <- do.call(rbindFill,
                                 lapply(split(spectra, f),
                                        .spectra_format_chromData))
              full_sp$chromIndex <- seq_len(nrow(full_sp))
              full_sp$chromSpectraIndex <- levels(f)
              rownames(full_sp) <- NULL
              object@spectra <- spectra
              callNextMethod(object, chromData = full_sp)
              }
          )

#' @rdname hidden_aliases
#' @importFrom methods callNextMethod
setMethod("show", "ChromBackendSpectra", function(object) {
    callNextMethod()
    cat("\nThe Spectra object contains", length(object@spectra), "spectra\n")
    if (object@inMemory) cat("\nPeaks data is cached in memory\n")
})

#' @rdname ChromBackendSpectra
chromSpectraIndex <- function(object) {
    if (!is(object, "ChromBackendSpectra"))
        stop("The object must be a 'ChromBackendSpectra' object.")
    chromData(object, columns = "chromSpectraIndex", drop = TRUE)
}

#' @rdname hidden_aliases
setMethod("factorize", "ChromBackendSpectra",
          function(object, factorize.by = c("msLevel", "dataOrigin"),...) {
              if (!all(factorize.by %in% chromVariables(object)))
                  stop("All 'factorize.by' variables must be in chromData.")
              if (!all(factorize.by %in% Spectra::spectraVariables(object@spectra)))
                  stop("All 'factorize.by' variables must be in the ",
                       "Spectra object.")
              object@chromData$chromSpectraIndex <- factor(do.call(
                  paste, c(object@chromData[, factorize.by], sep = "_")))
              object@spectra$chromSpectraIndex <- factor(
                  do.call(
                      paste,
                      c(as.list(Spectra::spectraData(object@spectra)[, factorize.by]),
                        sep = "_")),
                  levels = levels(chromSpectraIndex(object)))
              object
          })

#' @rdname hidden_aliases
#' @importMethodsFrom ProtGenerics backendParallelFactor
setMethod("backendParallelFactor", "ChromBackendSpectra", function(object, ...) {
    factor()
})

#' @rdname hidden_aliases
#' @export
setMethod("isReadOnly", "ChromBackendSpectra", function(object) TRUE)

#' @rdname hidden_aliases
setMethod("peaksData", "ChromBackendSpectra",
          function(object, columns = peaksVariables(object),
                   drop = FALSE, ...) {
              if (object@inMemory || !length(object)) return(callNextMethod())
              ## Ensure chromSpectraIndex only contains relevant levels needed
              valid_levels <- chromSpectraIndex(object)
              if (!all(levels(object@spectra$chromSpectraIndex) %in% valid_levels))
                  object@spectra$chromSpectraIndex <- factor(object@spectra$chromSpectraIndex,
                                                             levels = valid_levels)

              ## Process peaks data
              pd <- mapply(.process_peaks_data,
                           cd = split(chromData(object), valid_levels),
                           s = split(object@spectra, object@spectra$chromSpectraIndex),
                           MoreArgs = list(columns = columns,
                                           fun = object@summaryFun,
                                           drop = drop),
                           SIMPLIFY = FALSE)
              unlist(pd, use.names = FALSE, recursive = FALSE)
          })


#' @rdname hidden_aliases
setReplaceMethod("peaksData", "ChromBackendSpectra", function(object, value) {
    message("The `peaksData` slot will be modified but the changes will not",
            " affect the Spectra object.")
    object <- callNextMethod()
    object@inMemory <- TRUE
    object
})

#' @rdname hidden_aliases
setReplaceMethod("chromData", "ChromBackendSpectra", function(object, value) {
    message("Please keep in mind the 'ChromBackendSpectra' backend is read-only.",
            " The chromData slot will be modified but the changes will not",
            " affect the Spectra object. You will need to run `factorize()` to",
            " update the 'chromSpectraIndex' column.")
    callNextMethod()
})

#' @rdname hidden_aliases
#' @export
setMethod("supportsSetBackend", "ChromBackendSpectra",
          function(object, ...) FALSE)

#' @rdname hidden_aliases
#' @importMethodsFrom S4Vectors [ [[
#' @export
setMethod("[", "ChromBackendSpectra", function(x, i, j, ...) {
    if (!length(i)) return (ChromBackendSpectra())
    callNextMethod()
})

