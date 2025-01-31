#' @include helpers.R
#' @include hidden_aliases.R
#' @include ChromBackend.R
NULL

#' @title Chromatographic data backend for reading mzML files.
#'
#' @aliases [,ChromBackendMzR-method
#'
#' @name ChromBackendMzR
#'
#' @description
#' The `ChromBackendMzR` inherits all slots and methods from the base
#' `ChromBackendMemory` backend, providing additional functionality for
#' reading chromatographic data from mzML files.
#'
#' Contrarily to the `ChromBackendMemory` backend, the `ChromBackendMzR` backend
#' should have the *dataOrigin* chromatographic variables populated with the
#' file path of the mzML file from which the chromatographic data was read.
#'
#' It should be noted that the `ChromBackendMzR` backend is read-only, and does
#' not support direct modification of chromatographic data. However, it does
#' support `peaksData` slot replacement, which will modify the `peaksData` slot but
#' not the local mzML files. This is the  signaled by the "inMemory" being set
#' to TRUE.
#'
#' Lastly, implementing functionalities with the `ChromBackendMzR` backend
#' should be simplified as much as possible and reuse the methods already
#' implemented for `ChromBackendMemory` when possible.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param files A character vector of file paths to mzML files.
#'
#' @param object A `ChromBackendMzR` object.
#'
#' @param ... Additional parameters to be passed.
#'
#' @author Philippine Louail
#'
#' @exportClass ChromBackendMzR
#
NULL

#' @noRd
setClass("ChromBackendMzR",
         contains = "ChromBackendMemory",
         slots = c(inMemory = "logical"),
         prototype = prototype(chromData = fillCoreChromVariables(data.frame()),
                               peaksData = list(.EMPTY_PEAKS_DATA),
                               readonly = FALSE,
                               version = "0.1",
                               inMemory = FALSE))

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
setMethod("backendInitialize", "ChromBackendMzR",
          function(object, files = character(), BPPARAM = bpparam(),
                   ...) {
              if (!length(files)) return(object)
              if (!is.character(files))
                  stop("Parameter 'files' must be a character vector of file ",
                       "paths")
              files <- normalizePath(files, mustWork = FALSE)
              chromData <- do.call(rbindFill,
                                   bplapply(files,
                                            FUN = function(fl) {
                                                cbind(.mzR_format_chromData(fl))
                                                }, BPPARAM = BPPARAM))
              object <- callNextMethod(object, chromData = chromData, ...)
              object
          })

#' @rdname hidden_aliases
#' @importFrom utils head
#' @export
setMethod("show", "ChromBackendMzR", function(object) {
    callNextMethod()
    fls <- unique(object@chromData$dataOrigin)
    if (length(fls)) {
        to <- min(3, length(fls))
        cat("\nfile(s):\n", paste(basename(fls[seq_len(to)]), collapse = "\n"),
            "\n", sep = "")
        if (length(fls) > 3) cat(" ...", length(fls) - 3, "more files\n")
    }
    if (object@inMemory) cat("\nData is in memory in the peaksData slots\n")
})

#' @rdname hidden_aliases
#'
#' @importMethodsFrom ProtGenerics backendParallelFactor
setMethod("backendParallelFactor", "ChromBackendMzR", function(object, ...) {
    factor(dataOrigin(object), levels = unique(dataOrigin(object)))
    })


#' @rdname hidden_aliases
#' @export
setMethod("isReadOnly", "ChromBackendMzR", function(object) TRUE)

#' @rdname hidden_aliases
#' @importFrom BiocParallel bplapply
setMethod("peaksData",
          "ChromBackendMzR",
          function(object, columns = peaksVariables(object), drop = FALSE,
                   BPPARAM = SerialParam(), ...) {
              if (object@inMemory || !length(object)) return(callNextMethod())
              pv <- peaksVariables(object)
              if (!any(columns %in% pv))
                  stop("Some of the requested peaks variables are not",
                             " available")
              ret <- all(pv %in% columns)
              pd <- bplapply(split(object, f = dataOrigin(object)),
                             function(ob) {
                  chr <-.get_chrom_data(fl = unique(dataOrigin(ob)),
                                        idx = chromIndex(ob))
                  if (ret) chr
                  else lapply(chr, `[`, , columns, drop = drop)
              }, BPPARAM = BPPARAM)
              unlist(pd, recursive = FALSE, use.names = FALSE)
          })

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "ChromBackendMzR", function(object, value) {
    message("Please keep in mind the 'ChromBackendMzR' backend is read-only.",
            " The peaksData slot will be modified but the changes will not
            affect the local mzML files.")
    object <- callNextMethod()
    object@inMemory <- TRUE
    object
})

#' @rdname hidden_aliases
setReplaceMethod("chromData", "ChromBackendMzR", function(object, value) {
    message("Please keep in mind the 'ChromBackendMzR' backend is read-only.",
            " The chromData slot will be modified but the changes will not
            affect the local mzML files.")
    callNextMethod() # should there be extra check for storage and chromIndex ?
})

#' @rdname hidden_aliases
#' @export
setMethod("supportsSetBackend", "ChromBackendMzR", function(object, ...) FALSE)

#' @rdname hidden_aliases
setMethod("[", "ChromBackendMzR", function(x, i, j, ..., drop = FALSE) {
    if (!length(i)) return (ChromBackendMzR())
    callNextMethod()
})
