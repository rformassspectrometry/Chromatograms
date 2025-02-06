#' @include Chromatograms.R hidden_aliases.R ChromBackendMzR.R

#' @title Chromatographic peaks data
#'
#' @name peaksData
#'
#' @aliases peaksData
#' @aliases peaksVariables
#' @aliases filterPeaksData
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
#' @param ranges For `filterPeaksData()` : a `numeric` vector of paired values
#'        (upper and lower boundary) that define the ranges to filter the
#'        `object`. These paired values need to be in the same order as the
#'        `variables` parameter (see below).
#'
#' @param object A [Chromatograms] object.
#'
#' @param value For `rtime()` and `intensity()`: `numeric` vector with the
#'        values to replace the current values. The length of the vector must
#'        match the number of peaks data pairs in the `Chromatograms` object.
#'
#' @param variables For `filterPeaksData()`: `character` vector with the names
#'        of the peaks data variables to filter for. The list of available
#'        peaks data variables can be obtained with `peaksVariables()`.
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
          })

#' @rdname peaksData
#' @note read-only backend cannot be replaced in the chromatograms method BUT can
#' in the backend method.
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
    if (isReadOnly(object@backend))
        stop("Cannot replace peaks data in a read-only backend")
    rtime(object@backend) <- value
    object
})

#' @rdname peaksData
setMethod("intensity", signature = "Chromatograms", function(object, ...)
              peaksData(object, columns = "intensity", drop = TRUE))

#' @rdname peaksData
setReplaceMethod("intensity", signature = "Chromatograms", function(object, value) {
    if (isReadOnly(object@backend))
        stop("Cannot replace peaks data in a read-only backend")
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
                  object@processing, "Filter: remove peaks based ",
                  "on the varaibles: " , paste(variables, collapse = ", "),
                  "the ranges: ", paste(ranges, collapse = ", "),
                  "and the match condition: ", match)
              object
          })

