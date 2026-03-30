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
#' @aliases compareChromatograms
#' @aliases matchRtime
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
#'
#' @param FUN For `compareChromatograms()`: `function` to compute the
#'        similarity between two chromatograms from their aligned intensity
#'        vectors (as returned by `MAPFUN`). Must accept two numeric vectors
#'        as the first two arguments and return a single numeric value.
#'        Defaults to [stats::cor()] (Pearson correlation). Additional
#'        arguments can be passed via `...` (e.g., `method = "spearman"` for
#'        [stats::cor()]).
#'
#' @param MAPFUN For `compareChromatograms()`: `function` to align the
#'        retention times of two chromatograms before computing similarity.
#'        The function must accept two `data.frame`s (with columns `rtime`
#'        and `intensity`) and return a `list` with elements `x` and `y`
#'        (numeric vectors of equal length). Defaults to [matchRtime()],
#'        which linearly interpolates intensities onto a common retention
#'        time grid. Additional arguments can be passed via `...`.
#'
#' @param by For `compareChromatograms()`: `character` giving the name(s) of
#'        chromatogram variable(s) (columns in `chromData()`) to group
#'        chromatograms by before computing similarity. A named `list`
#'        of similarity matrices is returned, one per unique combination
#'        of the grouping variable(s). Defaults to `character()` (no
#'        grouping), which returns a single similarity matrix. Only used
#'        when `y` is missing.
#'
#' @param labels For `compareChromatograms()`: optional `character(1)` giving
#'        the name of a chromatogram variable (column in `chromData()`) whose
#'        values should be used as row and column names of the returned matrix.
#'        The column must contain unique values (within each group, if `by` is
#'        used). If `NULL` (the default), the matrix dimensions are unnamed.
#'        Only used when `y` is missing.
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
#' @param x For `lengths()` and `compareChromatograms()`: A [Chromatograms]
#'        object.
#'
#' @param y For `compareChromatograms()`: A [Chromatograms] object against
#'        which `x` is compared. If missing, each chromatogram in `x` is
#'        compared with each other chromatogram in `x`.
#'
#' @param SIMPLIFY For `compareChromatograms()`: `logical(1)` (default
#'        `TRUE`). If `TRUE` and either `x` or `y` has length 1, the result
#'        is simplified to a `numeric` vector instead of a `matrix`.
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
#' @section Compare Chromatograms:
#'
#' `compareChromatograms()` compares chromatograms following the same
#' two-step design as `compareSpectra()` in the *Spectra* package:
#'
#' 1. **Align** – `MAPFUN` (default [matchRtime()]) maps the two
#'    chromatograms onto a common retention-time grid and returns a
#'    `list(x, y)` of aligned intensity vectors.
#' 2. **Score** – `FUN` (default [stats::cor()], Pearson correlation)
#'    computes a single similarity value from the aligned intensities.
#'
#' If `y` is missing, each chromatogram in `x` is compared against every
#' other chromatogram in `x`; otherwise, each in `x` is compared with
#' each in `y`.
#'
#' The default `MAPFUN = matchRtime` linearly interpolates intensities
#' onto the union of both chromatograms' retention times within their
#' overlapping range. A custom `MAPFUN` can be supplied; it must accept
#' two `data.frame`s (with columns `rtime` and `intensity`) and return
#' a `list` with elements `x` and `y` (numeric vectors of equal length).
#'
#' The result is a `numeric` `matrix` with `length(x)` rows and
#' `length(y)` columns (or a symmetric `n x n` matrix for
#' self-comparison, with 1 on the diagonal). Pairs without overlapping
#' retention-time ranges, or with fewer than two aligned points, return
#' `NA`.
#'
#' If `SIMPLIFY = TRUE` (default) and `length(x)` or `length(y)` is one,
#' the result is simplified to a `numeric` vector.
#'
#' `matchRtime()` is the default retention-time alignment function used
#' by `compareChromatograms()`. It takes two `data.frame`s with columns
#' `rtime` and `intensity`, finds their overlapping retention-time range,
#' builds a common grid from the union of both chromatograms' time points
#' within that range, and linearly interpolates both chromatograms'
#' intensities onto the common grid using [stats::approx()]. The result
#' is a `list` with elements `x` and `y` (numeric vectors of equal
#' length). If the chromatograms do not overlap or either has fewer than
#' 2 data points, empty numeric vectors are returned.
#'
#' When `y` is missing, the `labels` parameter can be used to assign
#' meaningful row/column names to the output matrix from a `chromData()`
#' column (e.g., `"mz"` or a user-defined feature identifier). The column
#' must contain unique values (within each group, if `by` is used).
#'
#' When `y` is missing, the `by` parameter allows to split chromatograms
#' into groups based on one or more `chromData()` columns before computing
#' similarity. This is useful, for example, to compute separate similarity
#' matrices per `dataOrigin` or per `precursorMz`. When `by` is provided,
#' a named `list` of matrices is returned.
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
#' # Pairwise similarity of chromatogram intensity profiles
#' compareChromatograms(chr)
#' ## Use Spearman correlation (passed to cor() via ...)
#' compareChromatograms(chr, method = "spearman")
#'
#' # Use a chromData column as row/column labels
#' compareChromatograms(chr, labels = "mz")
#'
#' # Compare two Chromatograms objects
#' compareChromatograms(chr[1:2], chr[3])
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
matchRtime <- function(x, y, ...) {
    if (nrow(x) < 2L || nrow(y) < 2L)
        return(list(x = numeric(), y = numeric()))
    rt_min <- max(min(x$rtime), min(y$rtime))
    rt_max <- min(max(x$rtime), max(y$rtime))
    if (rt_min >= rt_max)
        return(list(x = numeric(), y = numeric()))
    common_rt <- sort(unique(c(
        x$rtime[x$rtime >= rt_min & x$rtime <= rt_max],
        y$rtime[y$rtime >= rt_min & y$rtime <= rt_max]
    )))
    if (length(common_rt) < 2L)
        return(list(x = numeric(), y = numeric()))
    list(
        x = approx(x$rtime, x$intensity, xout = common_rt, rule = 1)$y,
        y = approx(y$rtime, y$intensity, xout = common_rt, rule = 1)$y
    )
}

#' @rdname peaksData
#' @export
setMethod("compareChromatograms", 
          signature(x = "Chromatograms", y = "Chromatograms"),
    function(x, y, MAPFUN = matchRtime, FUN = cor, ..., SIMPLIFY = TRUE,
             BPPARAM = SerialParam()) {
        nx <- length(x)
        ny <- length(y)
        if (nx == 0L || ny == 0L)
            return(matrix(numeric(0), nx, ny))
        mat <- .compare_chromatograms(peaksData(x), peaksData(y),
                                      MAPFUN = MAPFUN, FUN = FUN,
                                      BPPARAM = BPPARAM, self = FALSE, ...)
        if (SIMPLIFY && (nx == 1L || ny == 1L))
            mat <- as.vector(mat)
        mat
    }
)

#' @rdname peaksData
#' @export
setMethod("compareChromatograms", signature(x = "Chromatograms", y = "missing"),
    function(x, y = NULL, MAPFUN = matchRtime, FUN = cor, ...,
             labels = NULL, by = character(), SIMPLIFY = TRUE,
             BPPARAM = SerialParam()) {
        if (length(by)) {
            if (!is.character(by))
                stop("'by' must be a character vector")
            missing_cols <- setdiff(by, chromVariables(x))
            if (length(missing_cols))
                stop("Column(s) '",
                     paste0(missing_cols, collapse = "', '"),
                     "' not found in chromData")
            grp <- interaction(chromData(x)[, by, drop = FALSE], drop = TRUE)
            lvls <- levels(grp)
            res <- vector("list", length(lvls))
            names(res) <- lvls
            for (g in lvls) {
                idx <- which(grp == g)
                sub_obj <- x[idx]
                labs <- .resolve_labels(sub_obj, labels)
                res[[g]] <- .compare_chromatograms(peaksData(sub_obj),
                    MAPFUN = MAPFUN, FUN = FUN, labels = labs,
                    BPPARAM = BPPARAM, self = TRUE, ...)
            }
            return(res)
        }
        n <- length(x)
        if (n == 0L) return(matrix(numeric(0), 0, 0))
        if (n == 1L)
            return(compareChromatograms(x, x, MAPFUN = MAPFUN,
                                        FUN = FUN, SIMPLIFY = FALSE,
                                        BPPARAM = BPPARAM, ...))
        labs <- .resolve_labels(x, labels)
        .compare_chromatograms(peaksData(x), MAPFUN = MAPFUN, FUN = FUN,
                               labels = labs, BPPARAM = BPPARAM, self = TRUE, ...)
    }
)

#' @rdname peaksData
#' @export
setMethod("peakBoundary", signature = "Chromatograms",
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
