#' @include Chromatograms.R hidden_aliases.R ChromBackendMzR.R

#' @title Chromatographic peaks data
#'
#' @name peaksData
#'
#' @aliases peaksData
#' @aliases peaksVariables
#' @aliases filterPeaksData
#' @aliases imputePeaksData
#'
#' @description
#'
#' As explained in the [Chromatograms] class documentation, the `Chromatograms`
#' object is a container for chromatographic data that includes chromatographic
#' peaks data (*retention time* and related intensity values, also referred to
#' as *peaks data variables* in the context of `Chromatograms`) and metadata of
#' individual chromatograms (so called *chromatograms variables*).
#'
#' The *peaks data variables* information can be accessed using the
#' `peaksData()` function. It is also possible to access specific peaks
#' variables using `$`.
#'
#' The peaks data can be accessed, replaced but also filtered/subsetted. Refer
#' to the sections below for more details.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param columns For `peaksData()`: optional `character` with column
#'        names (peaks variables) that should be included in the
#'        returned `list` of `data.frame`. By default, all columns are returned.
#'        Available variables can be found by calling `peaksVariables()` on the
#'        object.
#'
#' @param drop `logical(1)` For `peaksData()`, default to `FALSE`. If `TRUE`,
#'        and one column is called by the user, the method returns a list of
#'        vector of the single column requested.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param keep For `filterPeaksData()`: `logical(1)` defining whether to
#'        keep (`keep = TRUE`) or remove (`keep = FALSE`) the chromatographic
#'        peaks data that match the condition.
#'
#' @param match For `filterPeaksData()` : `character(1) ` defining whether the
#'        condition has to match for all provided `ranges` (`match = "all"`;
#'        the default), or for any of them (`match = "any"`).
#'
#' @param method For `imputePeaksData()`: `character(1)`: Imputation
#'        method ("linear", "spline", "gaussian", "loess")
#'
#' @param object A [Chromatograms] object.
#'
#' @param ranges For `filterPeaksData()` : a `numeric` vector of paired values
#'        (upper and lower boundary) that define the ranges to filter the
#'        `object`. These paired values need to be in the same order as the
#'        `variables` parameter (see below).
#'
#' @param value For `rtime()` and `intensity()`: `numeric` vector with the
#'        values to replace the current values. The length of the vector must
#'        match the number of peaks data pairs in the `Chromatograms` object.
#'
#' @param variables For `filterPeaksData()`: `character` vector with the names
#'        of the peaks data variables to filter for. The list of available
#'        peaks data variables can be obtained with `peaksVariables()`.
#'
#' @param sd For `imputePeaksData`: `numeric(1)`, for the gaussian method:
#'        Standard deviation for Gaussian kernel
#'        (only used if method == "gaussian")
#'
#' @param span For `imputePeaksData`: `numeric(1)`, for the loess method:
#'        Smoothing parameter (only used if method == "loess")
#'
#' @param window For `imputePeaksData`: `integer`, for the gaussian method:
#'        Half-width of Gaussian kernel window (e.g., 2 gives window size 5)
#'
#' @param ... Additional arguments passed to the method.
#'
#'
#' @section Filter Peaks Variables:
#'
#' Functions that filter a `Chromatograms`'s peaks data (i.e., `peaksData`).
#' These functions remove peaks data that do not meet the
#' specified conditions. If a chromatogram in a `Chromatograms` object is
#' filtered, only the corresponding peaks variable pairs (i.e., rows) in the
#' `peaksData` are removed, while the chromatogram itself remains in the object.
#'
#' The available functions to filter chromatographic peaks data include:
#'
#' - `filterPeaksData()`: Filters numerical peaks data variables based on the
#'   specified numerical `ranges` parameter. This method returns the same input
#'   `Chromatograms` object, but the filtering step is added to the processing
#'   queue. The filtered data will be reflected when the user accesses
#'   `peaksData`. This function does *not* reduce the number of chromatograms
#'   in the object, but it removes the specified peaks data (e.g., "rtime" and
#'   "intensity" pairs) from the `peaksData`.
#'
#' In the case of a read-only backend, (such as the [ChromBackendMzR]), the
#' replacement of the peaks data is not possible. The peaks data can be
#' filtered, but the filtered data will not be saved in the backend. This means
#' the original mzml files will not be affected by computations performed on
#' the [Chromatograms].
#'
#' @section Impute Peaks Variables:
#'
#' `imputePeaksData` will impute missing values in a `Chromatograms`'s peaks data
#' (i.e., `peaksData`). This functions replace missing peaks data values with
#' specified imputation methods using various methods such as linear
#' interpolation, spline interpolation, Gaussian kernel smoothing, or LOESS
#' smoothing. This method modifies the peaks data in place and returns the
#' same `Chromatograms` object with imputed values.
#'
#'
#'
#' @seealso [Chromatograms] for a general description of the `Chromatograms`
#'          object, and [chromData] for accessing,substituting and filtering
#'          chromatographic variables. For more information on the queuing
#'          of processings and parallelization for larger dataset processing
#'          see [processingQueue].
#'
#' @md
#'
#' @importFrom BiocParallel SerialParam bpparam
#'
#' @author Philippine Louail
#'
#' @return Refer to the individual function description for information on the
#' return value.
#'
#' @examples
#'
#' # Create a Chromatograms object
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     dataOrigin = c("mem1", "mem2", "mem3")
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
#' # Access peaks data
#' peaksData(chr)
#'
#' # Access specific peaks data variables
#' peaksData(chr, columns = "rtime")
#' rtime(chr)
#'
#' # Replace peaks data
#' rtime(chr)[[1]] <- c(1, 2, 3, 4)
#'
#' # Filter peaks data
#' filterPeaksData(chr, variables = "rtime", ranges = c(12.5, 13.5))
#'
NULL

#' @rdname peaksData
setMethod("imputePeaksData",
    signature = "Chromatograms",
    function(object,
             method = c("linear", "spline", "gaussian", "loess"),
             span = 0.3,
             sd = 1,
             window = 2,
             ...) {
        method <- match.arg(method)
        object <- addProcessing(object, imputePeaksData,
            method = method, span = span, sd = sd, window = window
        )
        object@processing <- .logging(
            .processing(object), "Impute: replace missing peaks data ",
            "using the '", method, "' method"
        )
        object
    }
)

#' @rdname peaksData
setMethod("filterPeaksData",
    signature = "Chromatograms",
    function(object, variables = character(),
    ranges = numeric(), match = c("any", "all"),
    keep = TRUE) {
        object <- addProcessing(object, filterPeaksData,
            variables = variables, ranges = ranges,
            match = match, keep = keep
        )
        object@processing <- .logging(
            .processing(object), "Filter: remove peaks based ",
            "on the variables: ", paste(variables, collapse = ", "),
            "the ranges: ", paste(ranges, collapse = ", "),
            "and the match condition: ", match
        )
        object
    }
)

#' @rdname peaksData
setMethod("intensity", signature = "Chromatograms", function(object, ...) {
    peaksData(object, columns = "intensity", drop = TRUE)
})

#' @rdname peaksData
setReplaceMethod("intensity",
                 signature = "Chromatograms",
                 function(object, value) {
                     if (isReadOnly(.backend(object))) {
                         stop("Cannot replace peaks data in a read-only backend")
                     }
                     intensity(object@backend) <- value
                     object
                 }
)

#' @rdname peaksData
#' @importFrom MsCoreUtils between
setMethod("peaksData",
          signature = "Chromatograms",
          function(object,
                   columns = peaksVariables(object),
                   f = processingChunkFactor(object),
                   BPPARAM = bpparam(), drop = FALSE, ...) {
              queue <- .processingQueue(object)
              if (length(queue)) {
                  bd <- .run_process_queue(.backend(object),
                                           queue = queue,
                                           f = f,
                                           BPPARAM = BPPARAM
                  )
                  return(peaksData(bd, columns = columns, drop = drop))
              }
              peaksData(.backend(object), columns = columns, drop = drop)
          }
)

#' @rdname peaksData
setReplaceMethod("peaksData",
                 signature = "Chromatograms",
                 function(object, value) {
                     if (isReadOnly(.backend(object))) {
                         stop("Cannot replace peaks data in a read-only backend")
                     }
                     peaksData(object@backend) <- value
                     object
                 }
)

#' @rdname peaksData
setMethod("peaksVariables", signature = "Chromatograms", function(object, ...) {
    peaksVariables(.backend(object))
})

#' @rdname peaksData
setMethod("rtime", signature = "Chromatograms", function(object, ...) {
    peaksData(object, columns = "rtime", drop = TRUE)
})

#' @rdname peaksData
setReplaceMethod("rtime",
                 signature = "Chromatograms",
                 function(object, value) {
                     if (isReadOnly(.backend(object))) {
                         stop("Cannot replace peaks data in a read-only backend")
                     }
                     rtime(object@backend) <- value
                     object
                 }
)

#' @rdname peaksData
setMethod("lengths", signature = "Chromatograms", function(x) {
    queue <- .processingQueue(object)
    if (length(queue)) {
        bd <- .run_process_queue(.backend(object),
                                 queue = queue,
                                 f = f,
                                 BPPARAM = BPPARAM
        )
        return(lengths(bd))
    }
    lengths(.backend(object))
})

