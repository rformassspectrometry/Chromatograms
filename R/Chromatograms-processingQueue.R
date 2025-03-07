#' @title Efficiently processing `Chromatograms` objects.
#'
#' @name processingQueue
#'
#' @description
#' The `processingQueue` of a `Chromatograms` object is a list of processing
#' steps (i.e., functions) that are stored within the object and applied only
#' when needed. This design allows data to be processed in a single step,
#' which is particularly useful for larger datasets. The processing queue
#' enables functions to be applied in a chunk-wise manner, facilitating
#' parallel processing and reducing memory demand.
#'
#' Since the peaks data can be quite large, a processing queue is used to
#' ensure efficiency. Generally, the processing queue is applied either
#' temporarily when calling `peaksData()` or permanently when calling
#' `applyProcessing()`. As explained below the processing efficiency can be
#' further improved by enabling chunk-wise processing.
#'
#' @section Apply Processing:
#'
#' The `applyProcessing()` function applies the processing queue to the backend
#' and returns the updated `Chromatograms` object. The processing queue is a
#' list of processing steps applied to the chromatograms data. Each element in
#' the list is a function that processes the chromatograms data. To apply
#' processing to the peaks data, the backend must be set to a non-read-only
#' backend using the `setBackend()` function.
#'
#' @section Parallel and Chunk-wise Processing of `Chromatograms`:
#'
#' Many operations on `Chromatograms` objects, especially those involving the
#' actual peaks data (see [peaksData]), support chunk-wise processing. This
#' involves splitting the `Chromatograms` into smaller parts (chunks) that are
#' processed iteratively. This enables parallel processing by data chunk and
#' reduces memory demand since only the peak data of the currently processed
#' subset is loaded into memory. Chunk-wise processing, which is disabled by
#' default, can be enabled by setting the processing chunk size of a
#' `Chromatograms` object using the `processingChunkSize()` function to a value
#' smaller than the length of the `Chromatograms` object. For example, setting
#' `processingChunkSize(chr) <- 1000` will cause any data manipulation
#' operation on `chr`, such as `filterPeaksData()`, to be performed in parallel
#' for sets of 1000 chromatograms in each iteration.
#'
#' Chunk-wise processing is particularly useful for `Chromatograms` objects
#' using an *on-disk* backend or for very large experiments. For small datasets
#' or `Chromatograms` using an in-memory backend, direct processing might be
#' more efficient. Setting the chunk size to `Inf` will disable chunk-wise
#' processing.
#'
#' Some backends may prefer a specific type of splitting and chunk-wise
#' processing. For example, the `ChromBackendMzR` backend needs to load MS data
#' from the original (mzML) files, so chunk-wise processing on a per-file basis
#' is ideal. The [Chromatograms::backendParallelFactor()] function for
#' `ChromBackend` allows backends to suggest a preferred data chunking by
#' returning a `factor` defining the respective data chunks. The
#' `ChromBackendMzR` returns a `factor` based on the *dataOrigin*
#' chromatograms variable. A `factor` of length 0 is returned if no particular
#' preferred splitting is needed. The suggested chunk definition will be used
#' if no finite `processingChunkSize()` is defined. Setting the
#' `processingChunkSize` overrides `backendParallelFactor`.
#'
#' Functions to configure parallel or chunk-wise processing:
#' - `processingChunkSize()`: Gets or sets the size of the chunks for parallel
#'   or chunk-wise processing of a `Chromatograms` object. With a value of
#'   `Inf` (the default), no chunk-wise processing will be performed.
#' - `processingChunkFactor()`: Returns a `factor` defining the chunks into
#'   which a `Chromatograms` object will be split for chunk-wise (parallel)
#'   processing. A `factor` of length 0 indicates that no chunk-wise processing
#'   will be performed.
#'
#' @note
#' Some backends might not support parallel processing. For these, the
#' `backendBpparam()` function will always return a `SerialParam()` regardless
#' of how parallel processing was defined.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param chunkSize `integer(1)` for `processingChunkFactor` defining the chunk
#'        size. The default is the value stored in the `Chromatograms` object's
#'        `processingChunkSize` slot.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param FUN For `addProcessing()`, a function to be added to the
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
#'         returns a `factor` defining the chunks into which `object` will be
#'         split for (parallel) chunk-wise processing or a `factor` of length 0
#'         if no splitting is defined.
#'
#' @author Johannes Rainer, Philippine Louail
#'
#' @importMethodsFrom ProtGenerics processingChunkSize processingChunkSize<-
#'   processingChunkFactor
#'
#' @examples
#' # Create a Chromatograms object
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     chromIndex = c(1L, 2L, 3L)
#' )
#'
#' pdata <- list(
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     ),
#'     data.frame(
#'         rtime = c(45.1, 46.2),
#'         intensity = c(100, 80.1)
#'     ),
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     )
#' )
#'
#' be <- backendInitialize(new("ChromBackendMemory"),
#'     chromData = cdata,
#'     peaksData = pdata
#' )
#'
#' chr <- Chromatograms(be)
#'
#' divide_intensities <- function(x, y, ...) {
#'     intensity(x) <- lapply(intensity(x), `/`, y)
#'     x
#' }
#'
#' ## Add the function to the procesing queue
#' chr <- addProcessing(chr, divide_intensities, y = 2)
#' chr
#'
#' # Apply the processing queue
#' chr <- applyProcessing(chr)
#'
NULL

#' @rdname processingQueue
#' @importMethodsFrom ProtGenerics applyProcessing
#' @export
setMethod(
    "applyProcessing", "Chromatograms",
    function(object, f = processingChunkFactor(object),
    BPPARAM = bpparam(), ...) {
        if (isReadOnly(object@backend)) {
            stop("Cannot apply processing to a read-only backend")
        }
        queue <- object@processingQueue
        if (!length(queue)) {
            return(object)
        }
        object@backend <- .run_process_queue(object@backend,
            queue = queue,
            f = f, BPPARAM, ...
        )
        object@processing <- .logging(
            object@processing,
            "Applied processing queue with ",
            length(queue),
            " steps"
        )
        object@processingQueue <- list()
        object
    }
)

#' @importFrom ProtGenerics ProcessingStep
#' @importMethodsFrom ProtGenerics addProcessing
#' @importClassesFrom ProtGenerics ProcessingStep
#' @rdname processingQueue
setMethod("addProcessing", "Chromatograms", function(object, FUN, ...) {
    if (missing(FUN)) {
        return(object)
    }
    object@processingQueue <- c(
        object@processingQueue,
        list(ProcessingStep(FUN, ARGS = list(...)))
    )
    validObject(object)
    object
})

#' @rdname processingQueue
#' @export
setMethod("processingChunkSize",
    signature = "Chromatograms",
    function(object, ...) object@processingChunkSize
)

#' @rdname processingQueue
#' @export
setReplaceMethod("processingChunkSize",
    signature = "Chromatograms",
    function(object, value) {
        object@processingChunkSize <- value
        object
    }
)

#' @rdname processingQueue
#' @importMethodsFrom ProtGenerics backendParallelFactor
#' @export
setMethod("processingChunkFactor",
    signature = "Chromatograms",
    function(object, chunkSize = processingChunkSize(object), ...) {
        lx <- length(object)
        if (is.finite(chunkSize) && chunkSize < lx) {
            fres <- as.factor(rep(seq_len(ceiling(lx / chunkSize)),
                each = chunkSize
            )[seq_len(lx)])
            return(fres)
        } else {
            backendParallelFactor(object@backend)
        }
    }
)
