#' @title Queuing processing `Chromatograms` objects for parallel processing
#'
#' @name processingQueue
#'
#' @description
#' The `processingQueue` of a Chromatograms is a list of processing steps that
#' are stored within the object and only applied when needed.
#' This was created so that the data can be processed in a single step and is very useful for larger datasets.
#' This is even more true as this processing queue will call function that can
#' be applied on the data in a chunk-wise manner. This allows for parallel processing
#' of the data and reduces the memory demand.
#'
#' @section Apply Processing:
#'
#' The `applyProcessing()` function applies the processing queue to the backend
#' and returns the updated `Chromatograms` object. The processing queue is a
#' list of processing steps that are applied to the chromatograms data. Each
#' element in the list is a function that processes the chromatograms data.
#' The user will need to set the backend to a non-read-only backend to apply
#' processing to the peaks data. This can be done by using the `setBackend()`
#' function.
#'
#' @section Parallel and chunk-wise processing of `Chromatograms`:
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
#' `ChromBackendMzR` returns for example a `factor` based on the *dataOrigin*
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
#'
#' @note
#'
#' Some backends might not support parallel processing at all.
#' For these, the `backendBpparam()` function will always return a
#' `SerialParam()` independently on how parallel processing was defined.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param chunkSize `integer(1)` for `processingChunkFactor` defining the chunk
#'        size. The defualt will be the value stored in the `Chromatograms`
#'        object's `processingChunkSize` slot.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param FUN For `addProcessing()` A function to be added to the
#'        `Chromatograms` object's processing queue.
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
NULL

#' @rdname processingQueue
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
#' @rdname processingQueue
setMethod("addProcessing", "Chromatograms", function(object, FUN, ...) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    validObject(object)
    object
})

#' @rdname processingQueue
#' @export
#'
setMethod("processingChunkSize", signature = "Chromatograms",
          function(object, ...) object@processingChunkSize)

#' @rdname processingQueue
#'
#' @export
setReplaceMethod("processingChunkSize", signature = "Chromatograms",
                 function(object, value) {
                     object@processingChunkSize <- value
                     object
                 })

#' @rdname processingQueue
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

