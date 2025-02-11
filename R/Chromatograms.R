#' @title The Chromatograms class to manage and access chromatographic data
#'
#' @name Chromatograms
#'
#' @aliases Chromatograms-class Chromatograms
#'
#' @description
#' The `Chromatograms` class encapsules chromatographic data and related
#' metadata. The chromatographic data is represented by a *backend* extending
#' the virtual [ChromBackend] class which provides the raw data to the
#' `Chromatograms` object. Different backends and their properties are
#' decribed in the [ChromBackend] class documentation.
#'
#' @section Creation of objects:
#'
#' `Chromatograms` objects can be created using the `Chromatograms()`
#' construction function.
#'
#' @section Data stored in a `Chromatograms` object:
#'
#' The `Chromatograms` object is a container for chromatographic data, which
#' includes peaks data (*retention time* and related intensity values, also
#' referred to as *peaks data variables* in the context of `Chromatograms`) and
#' metadata of individual chromatogram (so called *chromatograms variables*).
#' While a core set of chromatograms variables (the
#' `coreChromatogramsVariables()`) and peaks data variables (the
#' `corePeaksVariables()`) are guaranteed to be provided by a `Chromatograms`,
#' it is possible to add arbitrary variables to a `Chromatograms` object.
#'
#' The `Chromatograms` object is designed to contain chromatographic data of a
#' (large) set of chromatograms. The data is organized *linearly* and can be
#' thought of a list of chromatograms, i.e. each element in the `Chromatograms`
#' is one chromatogram.
#'
#' The *chromatograms variables* information in the `Chromatograms` object can
#' be accessed using the `chromData()` function. Specific chromatograms
#' variables can be accessed by either precising the `"columns"` parameter in
#' `chromData()` or using `$`. `chromData` can be accessed, replaced but
#' also filtered/subsetted. Refer to the [chromData] documentation for more
#' details.
#'
#' The *peaks data variables* information in the `Chromatograms` object can be
#' accessed using the `peaksData()` function. Specific peaks variables can be
#' accessed by either precising the `"columns"` parameter in `peaksData()` or
#' using `$`. `peaksData` can be accessed, replaced but also filtered/subsetted.
#' Refer to the [peaksData] documentation for more details.
#'
#' @section Processing of `Chromatograms` objects:
#'
#' Functions that process the chromatograms data in some ways can be applied to
#' the object either directly or by using the `processingQueue` mechanism. The
#' `processingQueue` is a list of processing steps that are stored within the
#' object and only applied when needed. This was created so that the data can be
#' processed in a single step and is very useful for larger datasets. This is
#' even more true as this processing queue will call function that can be applied
#' on the data in a chunk-wise manner. This allows for parallel processing of
#' the data and reduces the memory demand. To read more about the
#' `processingQueue`, and how to parallelize your processes, see the
#' [processingQueue] documentation.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param backend [ChromBackend] object providing the raw data for the
#'        `Chromatograms` object.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param name A `character` string specifying the name of the variable to
#'        access.
#'
#' @param object A [Chromatograms] object.
#'
#' @param processingQueue [list] a list of processing steps (i.e. functions) to
#'        be applied to the chromatographic data. The processing steps are
#'        applied in the order they are listed in the `processingQueue`.
#'
#' @param value The value to replace the variable with.
#'
#' @param x A [Chromatograms] object.
#'
#' @param ... Additional arguments.
#'
#' @md
#'
#' @exportClass Chromatograms
#'
#' @seealso [chromData] for a general description of the chromatographic
#'          metadata available in the object, as well as how to access, replace
#'          and subset them.
#'          [peaksData] for a general description of the chromatographic peaks
#'          data available in the object, as well as how to access, replace and
#'          subset them.
#'          [processingQueue] for more information on the queuing of
#'          processings and parallelization for larger dataset.
#'
#' @examples
#' ## Create a Chromatograms object
#' chroms <- Chromatograms(backend = ChromBackendMemory())
#'
NULL

#' The Chromatograms class
#'
#' The `Chromatograms` class is a container for chromatographic data.
#'
#' @slot backend [ChromBackend] the *backend* object providing the raw data for
#'       the `Chromatograms` object.
#'
#' @slot processingQueue [list] a list of processing steps to be applied to the
#'       `Chromatograms`. Each element in the list is a function that
#'       processes the data. The processing steps are applied in
#'       the order they are listed in the `processingQueue`.
#'
#' @slot processing [character] a character vector with the names of the
#'       processing steps that have been applied to the `Chromatograms` object.
#'       This is mainly used in the "show" method to display the processing
#'       steps that have been applied to the `Chromatograms` object.
#'
#' @slot processingChunkSize [numeric(1)] the number of chromatograms to be
#'       processed in a single chunk. This is useful for processing large
#'       data sets in smaller chunks to avoid memory issues.
#'
#' @slot version [character(1)] the version of the `Chromatograms` object.
#'
#' @noRd
setClass("Chromatograms",
         slots = c(
             backend = "ChromBackend",
             processingQueue = "list",
             processing = "character",
             processingChunkSize = "numeric",
             version = "character"),
         prototype = prototype(version = "0.1",
                               processingChunkSize = Inf)
         )

setValidity("Chromatograms", function(object) {
    msg <- character()
    if (!is(object@backend, "ChromBackend"))
        msg <- ("backend must be a ChromBackend object")
    if (!is.numeric(object@processingChunkSize) ||
        length(object@processingChunkSize) != 1)
        msg <- c(msg, "processingChunkSize must be a numeric value")
    msg <- c(msg, .valid_processing_queue(object@processingQueue))
    if (length(msg)) msg
    else TRUE
})

#' @rdname Chromatograms
#' @importFrom methods new
#' @export
Chromatograms <- function(backend = ChromBackendMemory(),
                          processingQueue = list(), ...) {
    new("Chromatograms", backend = backend,
        processingQueue = processingQueue, ...)
}

#' @rdname hidden_aliases
#'
#' @param object A [Chromatograms] object.
#' @importMethodsFrom methods show
#' @importFrom utils capture.output
#'
#' @exportMethod show
setMethod("show", "Chromatograms",
          function(object) {
              cat("Chromatographic data (", class(object)[1L], ") with ",
                  length(object@backend), " chromatograms in a ",
                  class(object@backend), " backend:\n", sep = "")
              if (length(object@backend)) {
                  txt <- capture.output(show(object@backend))
                  cat(txt[-1], sep = "\n")
              }
              if (length(object@processingQueue))
                  cat("Lazy evaluation queue:", length(object@processingQueue),
                      "processing step(s)\n")
              lp <- length(object@processing)
              if (lp) {
                  lps <- object@processing
                  if (lp > 3) {
                      lps <- lps[1:3]
                      lps <- c(lps, paste0("...", lp - 3, " more processings. ",
                                           "Use 'processingLog' to list all."))
                  }
                  cat("Processing:\n", paste(lps, collapse="\n "), "\n")
              }
          })

#' @rdname Chromatograms
#'
#' @note
#' This needs to be discussed, if we want for example to be able to set a
#' a backend to `ChromBackendMzR` we need to implement backendInitialize()
#' better. = Support peaksData and chromData as arguments AND have a way to
#' write .mzml files (which we do not have for chromatographic data).
#'
#' @importMethodsFrom ProtGenerics setBackend
#'
#' @exportMethod setBackend
setMethod(
    "setBackend", c("Chromatograms", "ChromBackend"),
    function(object, backend, f = processingChunkFactor(object),
             BPPARAM = SerialParam(), ...) {
        backend_class <- class(object@backend)
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        BPPARAM <- backendBpparam(backend, BPPARAM)# what
        if (!supportsSetBackend(backend))
            stop(class(backend), " does not support 'setBackend'")
        if (!length(f) || length(levels(f)) == 1 || !length(object))
            bd_new <- backendInitialize(backend, peaksData = peaksData(object),
                                    chromData = chromData(object))
        else {
            bd_new <- bplapply(
                split(object@backend, f = f),
                function(z, ...) {
                    backendInitialize(backend,
                                      peaksData = peaksData(z),
                                      chromData = chromData(z),
                                      BPPARAM = SerialParam())
                }, ..., BPPARAM = BPPARAM)
            bd_new <- backendMerge(bd_new)
        }
        object@backend <- bd_new
        object@processing <- .logging(object@processing,
                                           "Switch backend from ",
                                           backend_class, " to ",
                                           class(object@backend))
        object
        })

#' @rdname Chromatograms
#' @export
setMethod("$", signature = "Chromatograms", function(x, name) {
    x@backend[[name]]
})

#' @rdname Chromatograms
#' @export
setReplaceMethod("$", signature = "Chromatograms", function(x, name, value) {
    x@backend[[name]] <- value
    x
})
