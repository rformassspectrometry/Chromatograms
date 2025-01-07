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
#' The `Chromatograms` object is designed to contain chromatographic data of a
#' (large) set of chromatograms. The data is organized *linearly* and can be
#' thought of a list of chromatograms, i.e. each element in the `Chromatograms`
#' is one chomatogram.
#'
#' @section Creation of objects:
#'
#' `Chromatograms` objects can be created using the `Chromatograms()`
#' construction function.
#'
#' @param backend [ChromBackend] object providing the raw data for the
#'       `Chromatograms` object.
#'
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
             # processingQueueVariables = "character", # wait do we keep this ???
             processing = "character",
             processingChunkSize = "numeric",
             version = "character"),
         prototype = prototype(version = "0.1",
                               processingChunkSize = Inf)
         )

setValidity("Chromatograms", function(object) {
    if (!is(object@backend, "ChromBackend"))
        return("backend must be a ChromBackend object")
    if (!is.numeric(object@processingChunkSize) || length(object@processingChunkSize) != 1)
        return("processingChunkSize must be a numeric value")
    #check the processing queue
    return(TRUE)
})

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
#' @importFrom methods new
#' @export
Chromatograms <- function(backend = ChromBackendMemory(),
                         processingQueue = list(), ...) {
    new("Chromatograms", backend = backend,
        processingQueue = processingQueue, ...)
}
