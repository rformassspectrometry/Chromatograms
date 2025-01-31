#######################
##
## Chromatogram class
##
########################

#' @title The Chromatograms class to manage and access chromatographic data
#'
#' @name Chromatograms
#'
#' @aliases Chromatograms-class Chromatograms
#'
#' @description
#' The `Chromatograms` class encapsules chromatographic data and related metadata.
#' The chromatographic data is represented by a *backend* extending the virtual
#' [ChromBackend] class which provides the raw data to the `Chromatograms`
#' object. Different backends and their properties are explained in the
#' [ChromBackend] class documentation.
#'
#' @section Data stored in a `Chromatograms` object:
#'
#' The `Chromatograms` object is a container for chromatographic data that includes
#' peaks data (*retention time* and related intensity values,
#' also referred to as *peaks data variables* in the context of `Chromatograms`)
#' and metadata of individual chromatogram (so called
#' *chromatograms variables*). While a core set of chromatograms variables (the
#' `coreChromatogramsVariables()`) and peaks data variables (the
#' `corePeaksVariables()`) are guaranteed to be provided by a `Chromatograms`,
#' it is possible to add arbitrary additional chromatograms variables to a
#' `Chromatograms` object.
#'
#' The *chromatograms variables* information can be accessed using the
#' `chromData()` function. it is also possible to access specific
#' chromatograms variables using `$`. `chromData` can be accessed, replaced but
#' also filtered/subsetted. Refer to the [chromData] documentation for more
#' details.
#'
#' The `Chromatograms` object is designed to contain chromatographic data of a
#' (large) set of chromatograms. The data is organized *linearly* and can be
#' thought of a list of chromatograms, i.e. each element in the `Chromatograms`
#' is one chromatogram.
#'
#' @section Creation of objects:
#'
#' `Chromatograms` objects can be created using the `Chromatograms()`
#' construction function.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param backend [ChromBackend] object providing the raw data for the
#'       `Chromatograms` object.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param object A [Chromatograms] object.

#' @param processingQueue [list] a list of processing steps to be applied to the
#'        chromatograms data. Each element in the list is a function that
#'        processes the chromatograms data. The processing steps are applied in
#'        the order they are listed in the `processingQueue`.
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
#'          [processingQueue] for more information on the queuing
#'          of precessing and parallelization for larger dataset processing
#'          see .
#'
#' @examples
#' ## create a Chromatograms object
#' chroms <- Chromatograms(backend = ChromBackendMemory())
#'
NULL

#' The Chromatograms class
#'
#' The `Chromatograms` class is a container for chromatogram data and metadata
#' for LC-MS experiemnt.
#'
#' @slot backend [ChromBackend] the backend object providing the raw data for
#'       the `Chromatograms` object.
#'
#' @slot processingQueue [list] a list of processing steps to be applied to the
#'       chromatograms data. Each element in the list is a function that
#'       processes the chromatograms data. The processing steps are applied in
#'       the order they are listed in the `processingQueue`.
#'
#' @slot processing [character] a character vector with the names of the
#'       processing steps that have been applied to the chromatograms data.
#'       This is mainly used in the "show" method to display the processing
#'       steps that have been applied to the chromatograms data to the user.
#'
#' @slot processingChunkSize [numeric(1)] the number of chromatograms to be
#'       processed in a single chunk. This is useful for processing large
#'       chromatograms data sets in smaller chunks to avoid memory issues.
#'
#' @slot version [character(1)] the version of the `Chromatograms` object.
#'
#' @noRd
setClass("Chromatograms",
         slots = c(
             backend = "ChromBackend",
             processingQueue = "list",
             processing = "character", # could we rename to processing history ?
             processingChunkSize = "numeric",
             version = "character"),
         prototype = prototype(version = "0.1",
                               processingChunkSize = Inf)
         )

setValidity("Chromatograms", function(object) {
    msg <- character()
    if (!is(object@backend, "ChromBackend"))
        msg <- ("backend must be a ChromBackend object")
    if (!is.numeric(object@processingChunkSize) || length(object@processingChunkSize) != 1)
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
#' @importMethodsFrom ProtGenerics setBackend
#'
#' @exportMethod setBackend
setMethod(
    "setBackend", c("Chromatograms", "ChromBackend"),
    function(object, backend, f = processingChunkFactor(object),
             BPPARAM = bpparam(), ...) {
        backend_class <- class(object@backend)
        BPPARAM <- backendBpparam(object@backend, BPPARAM)
        BPPARAM <- backendBpparam(backend, BPPARAM)# wut
        if (!supportsSetBackend(backend))
            stop(class(backend), " does not support 'setBackend'")
        if (!length(f) || length(levels(f)) == 1 || !length(object))
            bd_new <- backendInitialize(backend, peaksData = peaksData(object),
                                    chromData = chromData(object), ...)
        else {
            bd_new <- bplapply(
                split(object@backend, f = f),
                function(z, ...) {
                    backendInitialize(backend,
                                      peaksData = peaksData(z),
                                      chromData = chromData(z),
                                      ...,
                                      BPPARAM = SerialParam())
                }, ..., BPPARAM = BPPARAM)
            bd_new <- backendMerge(bd_new)
        }
        object@backend <- bd_new
        object@processingQueue <- .logging(object@processing,
                                           "Switch backend from ",
                                           backend_class, " to ",
                                           class(object@backend))
        })

