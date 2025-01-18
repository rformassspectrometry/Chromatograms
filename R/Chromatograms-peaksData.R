#' @include Chromatograms.R hidden_aliases.R

#' @title Chromatographic peaks data
#'
#' @name peaksData
#'
#' @aliases peaksData
#' @aliases peaksVariables
#' @aliases filterPeaksData
#' @aliases applyProcessing
#'
#' @description
#'
#' As explained in the [`Chromatograms`] class documentation, the `Chromatograms`
#' object is a container for chromatogram data that includes chromatographic
#' peaks data (*retention time* and related intensity values, also referred to
#' as *peaks data variables* in the context of `Chromatograms`) and metadata of
#' individual chromatograms (so called *chromatograms variables*).
#'
#' The *peaks data variables* information can be accessed using the
#' `peaksData()` function. It is also possible to access specific
#' peaks variables using `$`.
#'
#' `peaks` can be accessed, replaced but also filtered/subsetted. Refer to
#' the sections below for more details.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param columns For `peaksData()` accessor: optional `character` with column
#'        names (peaks variables) that should be included in the
#'        returned `list` of `data.frame`. By default, all columns are returned.
#'
#' @param drop `logical(1)` default to `FALSE`. If `TRUE`, and one column is
#'        called by the user, the method returns a list of  vector of the
#'        single column requested.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param FUN For `addProcessing()` A function to be added to the
#'        `Chromatograms` object's processing queue.
#'
#' @param keep For `filterPeaksData()`: `logical(1)`
#'        defining whether to keep (`keep = TRUE`) or remove (`keep = FALSE`)
#'        the chromatogram data that match the condition.
#'
#' @param match For `filterPeaksData()` : `character(1) `
#'        defining whether the condition has to match for all provided
#'        `ranges` (`match = "all"`; the default), or for any of them
#'        (`match = "any"`).
#'
#' @param ranges For `filterPeaksData()` : a `numeric`
#'        vector of paired values (upper and lower boundary) that define the
#'        ranges to filter the `object`. These paired values need to be in the
#'        same order as the `variables` parameter (see below).
#'
#' @param object A [Chromatograms] object.
#'
#' @param value For `rtime()` and `intensity()`: `numeric` vector with the
#'        values to replace the current values.
#'
#' @param variables For `filterPeaksData()`: `character` vector with the names
#'        of the chromatogram variables to filter for. The list of available
#'        chromatogram variables can be obtained with `peaksVariables()`.
#'
#' @param ... Additional arguments passed to the method.
#'
#'
#' @section Filter Peaks Variables:
#'
#' Functions that filter `Chromatograms` based on peaks variables
#' (i.e., `peaksData`). These functions remove peaks data that do not meet the
#' specified conditions. If a chromatogram is filtered, only the corresponding
#' peaks variable pairs (i.e., rows) in the `peaksData` are removed, while the
#' chromatogram itself remains in the object.
#'
#' Since peaks data can be quite large, a processing queue is used to ensure
#' efficiency. The backend data remains unchanged until the processing queue is
#' applied (e.g., by running `applyProcessing()`). When calling `peaksData()`,
#' the processing queue is applied to the output, but the backend data is not
#' replaced.
#'
#' The available functions to filter chromatogram data include:
#'
#' - `filterPeaksData()`: Filters numerical peaks data variables based on the
#'   specified numerical `ranges`. This method returns the same input
#'   `Chromatograms` object, but the filtering step is added to the processing
#'   queue. The filtered data will be reflected when the user accesses
#'   `peaksData`. This function does *not* reduce the number of chromatograms
#'   in the object, but it removes the specified peaks data (e.g., "rtime" and
#'   "intensity" pairs) from the `peaksData`.

#'
#' @seealso [Chromatograms] for a general description of the `Chromatograms`
#'          object, and [chromData] for accessing,substituting and filtering
#'          chromatographic variables.
#'
#' @md
#'
#' @importFrom BiocParallel SerialParam bpparam
#'
#' @author Philippine Louail
NULL

#' @rdname peaksData
setMethod("peaksData",
          signature = "Chromatograms",
          function(object,
                   columns = peaksVariables(object),
                   f = processingChunkFactor(object),
                   BPPARAM = bpparam(), drop = FALSE, ...) {
              queue <- object@processingQueue
              if (length(queue)) {
                  bd <- .run_process_queue(object@backend,
                                           queue = queue,
                                           f = f,
                                           BPPARAM = BPPARAM)
                  return(peaksData(bd, columns = columns, drop = drop))
                  }
              peaksData(object@backend, columns = columns, drop = drop)
          }
)

#' @rdname peaksData
setReplaceMethod("peaksData", signature = "Chromatograms",
                 function(object, value) {
    if (isReadOnly(object@backend))
        stop("Cannot replace peaks data in a read-only backend")
    peaksData(object@backend) <- value
    object
})

#' @rdname peaksData
setMethod("peaksVariables", signature = "Chromatograms", function(object, ...)
              peaksVariables(object@backend))

#' @rdname peaksData
setMethod("rtime", signature = "Chromatograms", function(object, ...)
              peaksData(object, columns = "rtime", drop = TRUE))

#' @rdname peaksData
setReplaceMethod("rtime", signature = "Chromatograms", function(object, value) {
    rtime(object@backend) <- value
    object
})

#' @rdname peaksData
setMethod("intensity", signature = "Chromatograms", function(object, ...)
              peaksData(object, columns = "intensity", drop = TRUE))

#' @rdname peaksData
setReplaceMethod("intensity", signature = "Chromatograms", function(object, value) {
    intensity(object@backend) <- value
    object
})

#' @rdname peaksData
setMethod("filterPeaksData",
          signature = "Chromatograms",
          function(object, variables = character(),
                   ranges = numeric(), match = c("any", "all"),
                   keep = TRUE) {
              object <- addProcessing(object, filterPeaksData,
                                      variables = variables, ranges = ranges,
                                      match = match, keep = keep)
              object@processing <- .logging(
                  object@processing, "Filter: remove peaks ") #improve details of message
              object
          })

#' @rdname peaksData
#'
#' @importMethodsFrom ProtGenerics applyProcessing
#'
#' @export
setMethod("applyProcessing", "Chromatograms",
          function(object, f = processingChunkFactor(object),
                   BPPARAM = bpparam(), ...) {
    if (isReadOnly(object@backend))
       stop("Cannot apply processing to a read-only backend")
    queue <- object@processingQueue
    if (!length(queue)) return (object)
    object@backend  <- .run_process_queue(object@backend, queue = queue,
                                          f = f, BPPARAM, ...)
    object@processing <- .logging(object@processing,
                                  "Applied processing queue with ",
                                  length(queue),
                                  " steps")
    object@processingQueue <- list()
    object
})

#' @importFrom ProtGenerics ProcessingStep
#'
#' @importMethodsFrom ProtGenerics addProcessing
#'
#' @importClassesFrom ProtGenerics ProcessingStep
#'
#' @rdname peaksData
setMethod("addProcessing", "Chromatograms", function(object, FUN, ...) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    validObject(object)
    object
})


#' @title Parallel and chunk-wise processing of `Chromatograms`
#'
#' @rdname processingChunkSize
#'
#' @aliases processingChunkSize processingChunkSize<- processingChunkFactor
#'
#' @note This documentation is mostly a placeholder and will be updated when
#' the chunkwise implementation is finalized.
#'
#' @description
#'
#' Many operations on `Chromatograms` objects, specifically those working with
#' the actual peaks data (see [peaksData]), allow a chunk-wise processing in which
#' the `Chromatograms` is split into smaller parts (chunks) that are
#' iteratively processed. This enables parallel processing of the data (by
#' data chunk) and also reduces the memory demand since only the peak data
#' of the currently processed subset is loaded into memory and processed.
#' This chunk-wise processing, which is by default disabled, can be enabled
#' by setting the processing chunk size of a `Chromatograms` with the
#' `processingChunkSize()` function to a value which is smaller than the
#' length of the `Chromatograms` object. Setting `processingChunkSize(chr) <- 1000`
#' will cause any data manipulation operation on the `chr`, such as
#' `filterPeaksData()`, to be performed eventually in parallel for
#' sets of 1000 chromatograms in each iteration.
#'
#' Such chunk-wise processing is specifically useful for `Chromatograms` objects
#' using an *on-disk* backend or for very large experiments. For small data
#' sets or `Chromatograms` using an in-memory backend, a direct processing might
#' however be more efficient. Setting the chunk size to `Inf` will disable
#' the chunk-wise processing.
#'
#' For some backends a certain type of splitting and chunk-wise processing
#' might be preferable. The `ChromBackendMzR` backend for example needs to load
#' the MS data from the original (mzML) files, hence chunk-wise processing
#' on a per-file basis would be ideal. The [backendParallelFactor()] function
#' for `ChromBackend` allows backends to suggest a preferred splitting of the
#' data by returning a `factor` defining the respective data chunks. The
#' `ChromBackendMzR` returns for example a `factor` based on the *dataStorage*
#' chromatograms variable. A `factor` of length 0 is returned if no particular
#' preferred splitting should be performed. The suggested chunk definition
#' will be used if no finite `processingChunkSize()` is defined. Setting
#' the `processingChunkSize` overrides `backendParallelFactor`.
#'
#' Functions to configure parallel or chunk-wise processing:
#'
#' - `processingChunkSize()`: allows to get or set the size of the chunks for
#'   parallel processing or chunk-wise processing of a `Chromatograms` in general.
#'   With a value of `Inf` (the default) no chunk-wise processing will be
#'   performed.
#'
#' - `processingChunkFactor()`: returns a `factor` defining the chunks into
#'   which a `Chromatograms` will be split for chunk-wise (parallel) processing.
#'   A `factor` of length 0 indicates that no chunk-wise processing will be
#'   performed.
#'
#' @note
#'
#' Some backends might not support parallel processing at all.
#' For these, the `backendBpparam()` function will always return a
#' `SerialParam()` independently on how parallel processing was defined.
#'
#' @param chunkSize `integer(1)` for `processingChunkFactor` defining the chunk
#'        size. The defualt will be the value stored in the `Chromatograms`
#'        object's `processingChunkSize` slot.
#'
#' @param object A `Chromatograms` object.
#'
#' @param value `integer(1)` defining the chunk size.
#'
#' @param ... Additional arguments passed to the methods.
#'
#' @return `processingChunkSize()` returns the currently defined processing
#'         chunk size (or `Inf` if it is not defined). `processingChunkFactor()`
#'         returns a `factor` defining the chunks into which `object` will be split
#'         for (parallel) chunk-wise processing or a `factor` of length 0 if
#'         no splitting is defined.
#'
#' @author Johannes Rainer, Philippine Louail
#'
#' @importMethodsFrom ProtGenerics processingChunkSize processingChunkSize<- processingChunkFactor
#'
#' @export
setMethod("processingChunkSize", signature = "Chromatograms",
          function(object, ...) object@processingChunkSize)

#' @rdname processingChunkSize
#'
#' @export
setReplaceMethod("processingChunkSize", signature = "Chromatograms",
                 function(object, value) {
                     object@processingChunkSize <- value
                     object
                 })

#' @rdname processingChunkSize
#'
#' @export
setMethod("processingChunkFactor", signature = "Chromatograms",
          function(object, chunkSize = processingChunkSize(object), ...) {
              lx <- length(object)
              if (is.finite(chunkSize) && chunkSize < lx) {
                  fres <- as.factor(rep(1:ceiling(lx / chunkSize),
                                        each = chunkSize)[seq_len(lx)])
                  return(fres)
              }
              else backendParallelFactor(object@backend)
          })

