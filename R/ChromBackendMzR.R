#' @include helpers.R
#' @include hidden_aliases.R
#' @include ChromBackend.R
NULL

#' @title Chromatographic Data Backend for Reading mzML Files
#'
#' @name ChromBackendMzR
#'
#' @description
#' The `ChromBackendMzR` inherits all slots and methods from the base
#' `ChromBackendMemory` backend, providing additional functionality for reading
#' chromatographic data from mzML files.
#'
#' Unlike the `ChromBackendMemory` backend, the `ChromBackendMzR` backend
#' should have the *dataOrigin* chromatographic variables populated with the
#' file path of the mzML file from which the chromatographic data was read.
#'
#' Note that the `ChromBackendMzR` backend is read-only and does not support
#' direct modification of chromatographic data. However, it does support
#' `peaksData` slot replacement, which will modify the `peaksData` slot but not
#' the local mzML files. This is indicated by the "inMemory" slot being set to
#' TRUE.
#'
#' Implementing functionalities with the `ChromBackendMzR` backend should be
#' simplified as much as possible and reuse the methods already implemented for
#' `ChromBackendMemory` when possible.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#' @param files A character vector of file paths to mzML files.
#' @param object A `ChromBackendMzR` object.
#' @param ... Additional parameters to be passed.
#'
#' @author Philippine Louail
#'
#' @exportClass ChromBackendMzR
#'
#' @return Refer to the individual function description for information on the
#'         return value.
#'
#' @examples
#' library(mzR)
#' library(msdata)
#'
#' ## Load an mzML file
#' MRM_file <- system.file("proteomics", "MRM-standmix-5.mzML.gz",
#'     package = "msdata"
#' )
#'
#' ## Initialize the ChromBackendMzR object
#' be_empty <- ChromBackendMzR()
#' be <- backendInitialize(be_empty, files = MRM_file, BPPARAM = SerialParam())
#'
NULL

#' @noRd
setClass("ChromBackendMzR",
    contains = "ChromBackendMemory",
    slots = c(inMemory = "logical"),
    prototype = prototype(
        chromData = fillCoreChromVariables(data.frame()),
        peaksData = list(.EMPTY_PEAKS_DATA),
        readonly = FALSE,
        version = "0.1",
        inMemory = FALSE
    )
)

#' @rdname ChromBackendMzR
#' @export ChromBackendMzR
ChromBackendMzR <- function() {
    .check_mzR_package()
    new("ChromBackendMzR")
}

#' @rdname ChromBackendMzR
#' @importFrom methods callNextMethod
#' @importFrom MsCoreUtils rbindFill
#' @export
setMethod(
    "backendInitialize", "ChromBackendMzR",
    function(object, files = character(), BPPARAM = bpparam(), ...) {
        if (!length(files)) {
            return(object)
        }
        if (!is.character(files)) {
            stop(
                "Parameter 'files' must be a character vector of ",
                "file paths"
            )
        }
        files <- normalizePath(files, mustWork = FALSE)
        chromData <- do.call(
            rbindFill,
            bplapply(files, FUN = function(fl) {
                cbind(.mzR_format_chromData(fl))
            }, BPPARAM = BPPARAM)
        )
        callNextMethod(object, chromData = chromData, ...)
    }
)

#' @rdname hidden_aliases
#' @export
setMethod("show", "ChromBackendMzR", function(object) {
    callNextMethod()
    fls <- unique(dataOrigin(object))
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n", paste(basename(fls[seq_len(to)]), collapse = "\n"),
            "\n",
            sep = ""
        )
        if (length(fls) > 3) cat(" ...", length(fls) - 3, "more files\n")
    }
    if (object@inMemory) cat("\nPeaks data is cached in memory\n")
})

#' @rdname hidden_aliases
#' @importMethodsFrom ProtGenerics backendParallelFactor
setMethod("backendParallelFactor", "ChromBackendMzR", function(object, ...) {
    factor(dataOrigin(object), levels = unique(dataOrigin(object)))
})

#' @rdname hidden_aliases
#' @export
setMethod("isReadOnly", "ChromBackendMzR", function(object) TRUE)

#' @rdname hidden_aliases
#' @importFrom BiocParallel bplapply
setMethod("peaksData", "ChromBackendMzR",
          function(object, columns = peaksVariables(object), drop = FALSE,
                   BPPARAM = SerialParam(), ...) {
              if (.inMemory(object) || !length(object)) return(callNextMethod())
              pv <- peaksVariables(object)
              if (!any(columns %in% pv))
                  stop("Some of the requested peaks variables are not",
                       " available")
              ret <- all(pv %in% columns)
              f <- factor(dataOrigin(object),
                          levels = unique(dataOrigin(object)))
              pd <- bplapply(split(object, f = f),
                 function(ob) {
                     chr <- .get_chrom_data(fl = .chromData(ob)$dataOrigin[1L],
                                            idx = chromIndex(ob))
                     if (ret) chr
                     else lapply(chr, `[`, , columns, drop = drop)
                 }, BPPARAM = BPPARAM)
              unsplit(pd, f = f)
          })

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "ChromBackendMzR", function(object, value) {
    message(
        "Please keep in mind the 'ChromBackendMzR' backend is read-only.",
        " The peaksData slot will be modified but the changes will not",
        " affect the local mzML files."
    )
    object <- callNextMethod()
    object@inMemory <- TRUE
    object
})

#' @rdname hidden_aliases
setReplaceMethod("chromData", "ChromBackendMzR", function(object, value) {
    message(
        "Please keep in mind the 'ChromBackendMzR' backend is read-only.",
        " The chromData slot will be modified but the changes will not",
        " affect the local mzML files."
    )
    callNextMethod()
})

#' @rdname hidden_aliases
#' @export
setMethod("supportsSetBackend", "ChromBackendMzR", function(object, ...) FALSE)

#' @rdname hidden_aliases
#' @importMethodsFrom S4Vectors [ [<-
setMethod("[", "ChromBackendMzR", function(x, i, j, ...) {
    if (!length(i)) {
        return(ChromBackendMzR())
    }
    callNextMethod()
})
