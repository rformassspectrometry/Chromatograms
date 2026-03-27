#' @include Chromatograms.R hidden_aliases.R ChromBackendMzR.R

#' @title Chromatographic peaks data
#'
#' @name peaksData
#'
#' @aliases peaksData
#' @aliases peaksVariables
#' @aliases filterPeaksData
#' @aliases imputePeaksData
#' @aliases peakBoundary
#' @aliases correlate
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
#'        method ("linear", "spline", "gaussian", "loess").
#'        For `correlate()`: `character(1)`: Correlation method
#'        ("pearson", "spearman", "kendall"). Default is "pearson".
#'
#' @param FUN For `correlate()`: optional custom similarity function. If
#'        provided, it must accept two numeric vectors (interpolated
#'        intensities) as the first two arguments and return a single numeric
#'        value. When `FUN` is provided, the `method` parameter is ignored.
#'        Additional arguments can be passed via `...`.
#'
#' @param by For `correlate()`: `character` giving the name(s) of
#'        chromatogram variable(s) (columns in `chromData()`) to group
#'        chromatograms by before computing correlation. A named `list`
#'        of correlation matrices is returned, one per unique combination
#'        of the grouping variable(s). Defaults to `character()` (no
#'        grouping), which returns a single correlation matrix.
#'
#' @param labels For `correlate()`: optional `character(1)` giving the name of
#'        a chromatogram variable (column in `chromData()`) whose values should
#'        be used as row and column names of the returned matrix. The column
#'        must contain unique values (within each group, if `by` is used).
#'        If `NULL` (the default), the matrix dimensions are unnamed.
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
#' @param extrapolate For `imputePeaksData`: `logical(1)` (default `FALSE`).
#'        If `TRUE`, missing values at the beginning and end of a chromatogram
#'        (outside the range of observed values) will be extrapolated. If
#'        `FALSE`, only interpolation is performed and leading/trailing `NA`
#'        values remain `NA`.
#'
#' @param sd For `imputePeaksData`: `numeric(1)`, for the gaussian method:
#'        Standard deviation for Gaussian kernel
#'        (only used if method == "gaussian")
#'
#' @param span For `imputePeaksData`: `numeric(1)`, for the loess method:
#'        Smoothing parameter (only used if method == "loess")
#'
#' @param threshold For `peakBoundary()`: `numeric(1)` (default `0.1`).
#'        Fraction of the peak height above baseline used as a fallback
#'        cut-off when valley-based boundaries are not suitable. Must be
#'        `>= 0` and `< 1`.
#'
#' @param baselineThreshold For `peakBoundary()`: `numeric(1)` (default
#'        `0.1`). Fraction of the peak height above the baseline. Valley
#'        positions returned by `MsCoreUtils::valleys()` are accepted only if
#'        the intensity at the valley is at or below
#'        `baseline + peak_height * baselineThreshold`. Must be `>= 0` and
#'        `< 1`.
#'
#' @param baselineQuantile For `peakBoundary()`: `numeric(1)` (default
#'        `0.1`). Quantile of the intensity distribution used as the
#'        baseline estimate. Must be `>= 0` and `<= 1`.
#'
#' @param window For `imputePeaksData`: `integer`, for the gaussian method:
#'        Half-width of Gaussian kernel window (e.g., 2 gives window size 5)
#'
#' @param x For `lengths()`: A `Chromatograms` object.
#'
#' @param ... Additional arguments passed to the method.
#'
#'
#' @section Filter Peaks Variables:
#'
#' Functions that filter a `Chromatograms`'s peaks data (i.e., `@peaksData`).
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
#' the original mzML files will not be affected by computations performed on
#' the [Chromatograms].
#'
#' @section Impute Peaks Variables:
#'
#' `imputePeaksData` will impute missing values in a `Chromatograms`'s peaks data
#' (i.e., `@peaksData`). This functions replace missing peaks data values with
#' specified imputation methods using various methods such as linear
#' interpolation, spline interpolation, Gaussian kernel smoothing, or LOESS
#' smoothing. This method modifies the peaks data in place and returns the
#' same `Chromatograms` object with imputed values.
#'
#'
#'
#' @section Peak Boundary Detection:
#'
#' `peakBoundary()` determines the retention time boundaries of the tallest
#' peak in each chromatogram. The function uses `MsCoreUtils::valleys()` to
#' locate the valleys (local minima) flanking the apex. If the valley
#' intensities exceed a baseline-relative threshold (controlled by
#' `baselineThreshold`), it falls back to a threshold-based boundary search
#' using `threshold`. The baseline is estimated as the `baselineQuantile`
#' quantile of the chromatogram's intensity values.
#' The result is a `matrix` with one row per
#' chromatogram and columns `left_boundary` and `right_boundary`
#' (retention times). Chromatograms that are empty, have fewer than 3 data
#' points, contain only `NA` or all-zero intensities return `NA` for both
#' boundaries.
#'
#'
#' @section Correlation:
#'
#' `correlate()` computes pairwise similarity (by default Pearson correlation)
#' between the intensity profiles of all chromatograms in the object.
#' Because chromatograms may have different retention time points, intensities
#' are linearly interpolated onto a common retention time grid (the union of
#' both chromatograms' retention times within their overlapping range) before
#' computing the similarity.
#'
#' The result is a square numeric `matrix` of dimensions `n x n`, where `n`
#' is the number of chromatograms. Diagonal elements are always `1`. Pairs
#' without overlapping retention time ranges, or with fewer than two data
#' points each, return `NA`.
#'
#' A custom similarity function can be passed via the `FUN` parameter. It
#' must accept two numeric vectors and return a single numeric value. When
#' `FUN` is provided, the `method` parameter is ignored.
#'
#' The `labels` parameter can be used to assign meaningful row/column names
#' to the output matrix from a `chromData()` column (e.g., `"mz"` or a
#' user-defined feature identifier). The column must contain unique values
#' (within each group, if `by` is used).
#'
#' The `by` parameter allows to split chromatograms into groups based on one
#' or more `chromData()` columns before computing correlation. This is useful,
#' for example, to compute separate correlation matrices per `dataOrigin` or
#' per `precursorMz`. When `by` is provided, a named `list` of matrices is
#' returned.
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
#' # Pairwise correlation of chromatogram intensity profiles
#' correlate(chr)
#' correlate(chr, method = "spearman")
#'
#' # Use a chromData column as row/column labels
#' correlate(chr, labels = "mz")
#'
NULL

#' @rdname peaksData
setMethod(
    "imputePeaksData", signature = "Chromatograms",
    function(object, method = c("linear", "spline", "gaussian", "loess"),
             span = 0.3, sd = 1, window = 2, extrapolate = FALSE, ...) {
        method <- match.arg(method)
        object <- addProcessing(object, imputePeaksData,
            method = method, span = span, sd = sd,
            window = window, extrapolate = extrapolate)
        object@processing <- .logging(
            .processing(object),
            "Impute: replace missing peaks data ",
            "using the '", method, "' method",
            if (extrapolate) " (with extrapolation)" else "")
        object
    }
)

#' @rdname peaksData
setMethod(
    "filterPeaksData", signature = "Chromatograms",
    function(object, variables = character(), ranges = numeric(),
             match = c("any", "all"), keep = TRUE) {
        object <- addProcessing(object, filterPeaksData,
            variables = variables, ranges = ranges,
            match = match, keep = keep)
        object@processing <- .logging(
            .processing(object),
            "Filter: remove peaks based ",
            "on the variables: ",
            paste(variables, collapse = ", "),
            "the ranges: ",
            paste(ranges, collapse = ", "),
            "and the match condition: ", match)
        object
    }
)

#' @rdname peaksData
setMethod("intensity", signature = "Chromatograms", function(object, ...) {
  if (!length(.processingQueue(object)))
    return(intensity(.backend(object)))
  peaksData(object, columns = "intensity", drop = TRUE)
})

#' @rdname peaksData
setReplaceMethod(
    "intensity", signature = "Chromatograms",
    function(object, value) {
        if (isReadOnly(.backend(object)))
            stop("Cannot replace peaks data in a read-only backend")
        intensity(object@backend) <- value
        object
    }
)

#' @rdname peaksData
#' @importFrom MsCoreUtils between
setMethod(
    "peaksData", signature = "Chromatograms",
    function(object, columns = peaksVariables(object),
             f = processingChunkFactor(object),
             BPPARAM = bpparam(), drop = FALSE, ...) {
        queue <- .processingQueue(object)
        if (length(queue)) {
            bd <- .run_process_queue(.backend(object),
                queue = queue, f = f, BPPARAM = BPPARAM)
            return(peaksData(bd, columns = columns, drop = drop))
        }
        peaksData(.backend(object), columns = columns, drop = drop)
    }
)

#' @rdname peaksData
setReplaceMethod(
    "peaksData", signature = "Chromatograms",
    function(object, value) {
        if (isReadOnly(.backend(object)))
            stop("Cannot replace peaks data in a read-only backend")
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
  if (!length(.processingQueue(object)))
    return(rtime(.backend(object)))
  peaksData(object, columns = "rtime", drop = TRUE)
})

#' @rdname peaksData
setReplaceMethod("rtime", signature = "Chromatograms", function(object, value) {
  if (isReadOnly(.backend(object))) {
    stop("Cannot replace peaks data in a read-only backend")
  }
  rtime(object@backend) <- value
  object
})

#' @rdname peaksData
setMethod("lengths", signature = "Chromatograms", function(x) {
    queue <- .processingQueue(x)
    if (length(queue)) {
        f <- processingChunkFactor(x)
        bd <- .run_process_queue(.backend(x),
            queue = queue, f = f, BPPARAM = bpparam())
        return(lengths(bd))
    }
    lengths(.backend(x))
})

#' @rdname peaksData
#' @importFrom stats cor
#' @export
setMethod(
    "correlate", signature = "Chromatograms",
    function(object, method = c("pearson", "spearman", "kendall"),
             FUN = NULL, labels = NULL, by = character(), ...) {
        method <- match.arg(method)
        if (length(by)) {
            if (!is.character(by))
                stop("'by' must be a character vector")
            missing_cols <- setdiff(by, chromVariables(object))
            if (length(missing_cols))
                stop("Column(s) '", paste0(missing_cols, collapse = "', '"),
                     "' not found in chromData")
            grp <- interaction(chromData(object)[, by, drop = FALSE],
            drop = TRUE)
            lvls <- levels(grp)
            res <- vector("list", length(lvls))
            names(res) <- lvls
            for (g in lvls) {
                idx <- which(grp == g)
                sub_obj <- object[idx]
                labs <- .resolve_labels(sub_obj, labels)
                res[[g]] <- .correlate_matrix(peaksData(sub_obj),
                    method = method, FUN = FUN, labels = labs, ...)
            }
            return(res)
        }
        n <- length(object)
        if (n == 0L) return(matrix(numeric(0), 0, 0))
        labs <- .resolve_labels(object, labels)
        .correlate_matrix(peaksData(object), method = method, FUN = FUN,
                          labels = labs, ...)
    }
)

#' @rdname peaksData
#' @export
setMethod(
    "peakBoundary", signature = "Chromatograms",
    function(object, threshold = 0.1, baselineThreshold = 0.1,
             baselineQuantile = 0.1, ...) {
        if (!is.numeric(threshold) || length(threshold) != 1L ||
            is.na(threshold))
            stop("'threshold' must be a single non-missing numeric value.")
        if (threshold < 0 || threshold >= 1)
            stop("'threshold' must be >= 0 and < 1.")
        if (!is.numeric(baselineThreshold) ||
            length(baselineThreshold) != 1L || is.na(baselineThreshold))
            stop("'baselineThreshold' must be a single non-missing ",
                 "numeric value.")
        if (baselineThreshold < 0 || baselineThreshold >= 1)
            stop("'baselineThreshold' must be >= 0 and < 1.")
        if (!is.numeric(baselineQuantile) ||
            length(baselineQuantile) != 1L || is.na(baselineQuantile))
            stop("'baselineQuantile' must be a single non-missing ",
                 "numeric value.")
        if (baselineQuantile < 0 || baselineQuantile > 1)
            stop("'baselineQuantile' must be >= 0 and <= 1.")
        pd <- peaksData(object)
        res <- vapply(pd, function(p)
            .peak_boundary_one(p$rtime, p$intensity,
                               threshold = threshold,
                               baselineThreshold = baselineThreshold,
                               baselineQuantile = baselineQuantile),
            numeric(2))
        res <- t(res)
        colnames(res) <- c("left_boundary", "right_boundary")
        res
    }
)
