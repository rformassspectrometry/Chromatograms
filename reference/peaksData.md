# Chromatographic peaks data

As explained in the
[Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
class documentation, the `Chromatograms` object is a container for
chromatographic data that includes chromatographic peaks data
(*retention time* and related intensity values, also referred to as
*peaks data variables* in the context of `Chromatograms`) and metadata
of individual chromatograms (so called *chromatograms variables*).

The *peaks data variables* information can be accessed using the
`peaksData()` function. It is also possible to access specific peaks
variables using `$`.

The peaks data can be accessed, replaced but also filtered/subsetted.
Refer to the sections below for more details.

## Usage

``` r
# S4 method for class 'Chromatograms'
imputePeaksData(
  object,
  method = c("linear", "spline", "gaussian", "loess"),
  span = 0.3,
  sd = 1,
  window = 2,
  extrapolate = FALSE,
  ...
)

# S4 method for class 'Chromatograms'
filterPeaksData(
  object,
  variables = character(),
  ranges = numeric(),
  match = c("any", "all"),
  keep = TRUE
)

# S4 method for class 'Chromatograms'
intensity(object, ...)

# S4 method for class 'Chromatograms'
intensity(object) <- value

# S4 method for class 'Chromatograms'
peaksData(
  object,
  columns = peaksVariables(object),
  f = processingChunkFactor(object),
  BPPARAM = bpparam(),
  drop = FALSE,
  ...
)

# S4 method for class 'Chromatograms'
peaksData(object) <- value

# S4 method for class 'Chromatograms'
peaksVariables(object, ...)

# S4 method for class 'Chromatograms'
rtime(object, ...)

# S4 method for class 'Chromatograms'
rtime(object) <- value

# S4 method for class 'Chromatograms'
lengths(x)

matchRtime(x, y, tolerance = Inf, ...)

# S4 method for class 'Chromatograms,Chromatograms'
compareChromatograms(
  x,
  y,
  MAPFUN = matchRtime,
  FUN = cor,
  ...,
  minPeaks = 4L,
  BPPARAM = SerialParam()
)

# S4 method for class 'Chromatograms,missing'
compareChromatograms(
  x,
  y = NULL,
  MAPFUN = matchRtime,
  FUN = cor,
  ...,
  minPeaks = 4L,
  labelsColumn = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'Chromatograms'
peakBoundary(
  object,
  threshold = 0.1,
  baselineThreshold = 0.1,
  baselineQuantile = 0.1,
  ...
)
```

## Arguments

- object:

  A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object.

- method:

  For `imputePeaksData()`: `character(1)`: Imputation method ("linear",
  "spline", "gaussian", "loess").

- span:

  For `imputePeaksData`: `numeric(1)`, for the loess method: Smoothing
  parameter (only used if method == "loess")

- sd:

  For `imputePeaksData`: `numeric(1)`, for the gaussian method: Standard
  deviation for Gaussian kernel (only used if method == "gaussian")

- window:

  For `imputePeaksData`: `integer`, for the gaussian method: Half-width
  of Gaussian kernel window (e.g., 2 gives window size 5)

- extrapolate:

  For `imputePeaksData`: `logical(1)` (default `FALSE`). If `TRUE`,
  missing values at the beginning and end of a chromatogram (outside the
  range of observed values) will be extrapolated. If `FALSE`, only
  interpolation is performed and leading/trailing `NA` values remain
  `NA`.

- ...:

  Additional arguments passed to the method.

- variables:

  For `filterPeaksData()`: `character` vector with the names of the
  peaks data variables to filter for. The list of available peaks data
  variables can be obtained with `peaksVariables()`.

- ranges:

  For `filterPeaksData()` : a `numeric` vector of paired values (upper
  and lower boundary) that define the ranges to filter the `object`.
  These paired values need to be in the same order as the `variables`
  parameter (see below).

- match:

  For `filterPeaksData()` : `character(1) ` defining whether the
  condition has to match for all provided `ranges` (`match = "all"`; the
  default), or for any of them (`match = "any"`).

- keep:

  For `filterPeaksData()`: `logical(1)` defining whether to keep
  (`keep = TRUE`) or remove (`keep = FALSE`) the chromatographic peaks
  data that match the condition.

- value:

  For
  [`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  and
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `numeric` vector with the values to replace the current values. The
  length of the vector must match the number of peaks data pairs in the
  `Chromatograms` object.

- columns:

  For `peaksData()`: optional `character` with column names (peaks
  variables) that should be included in the returned `list` of
  `data.frame`. By default, all columns are returned. Available
  variables can be found by calling `peaksVariables()` on the object.

- f:

  `factor` defining the grouping to split the `Chromatograms` object.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- drop:

  `logical(1)` For `peaksData()`, default to `FALSE`. If `TRUE`, and one
  column is called by the user, the method returns a list of vector of
  the single column requested.

- x:

  For [`lengths()`](https://rdrr.io/r/base/lengths.html) and
  `compareChromatograms()`: A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object. For `matchRtime()`: a `data.frame` with columns `rtime` and
  `intensity` representing the first chromatogram.

- y:

  For `compareChromatograms()`: A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object against which `x` is compared. If missing, each chromatogram in
  `x` is compared with each other chromatogram in `x`. For
  `matchRtime()`: a `data.frame` with columns `rtime` and `intensity`
  representing the second chromatogram.

- tolerance:

  For `matchRtime()`: `numeric(1)` (default `Inf`). Maximum RT
  difference between two measured points to be considered a match.
  Controls both the overlap detection and the shared RT grid. Lower
  values prevent a peak from being compared against a long interpolated
  gap in the other chromatogram. Use `Inf` (the default) to consider all
  RT points as matching. Can be forwarded via `...` in
  `compareChromatograms()`.

- MAPFUN:

  For `compareChromatograms()`: `function` to align the retention times
  of two chromatograms before computing similarity. Must accept two
  `data.frame`s (with columns `rtime` and `intensity`) and return a
  `list` with elements `x` and `y`: numeric vectors of equal length
  containing the aligned intensities of the first and second
  chromatogram respectively, interpolated onto a common retention-time
  grid. Defaults to `matchRtime()`. Additional arguments can be passed
  via `...`.

- FUN:

  For `compareChromatograms()`: `function` to compute the similarity
  between two chromatograms from their aligned intensity vectors (as
  returned by `MAPFUN`). Must accept two numeric vectors as the first
  two arguments and return a single numeric value. Defaults to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) (Pearson
  correlation). Additional arguments can be passed via `...` (e.g.,
  `method = "spearman"` for
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)).

- minPeaks:

  For `compareChromatograms()`: `integer(1)` (default `4L`). Minimum
  number of overlapping retention-time points (as returned by `MAPFUN`)
  required to compute a similarity score. Pairs whose retention-time
  overlap contains fewer than `minPeaks` points return `NA` in the score
  layer; the actual overlap count is still recorded in the `n_peaks`
  layer. Setting `minPeaks = 2L` recovers the behaviour of always
  computing a score whenever at least two points overlap.

- labelsColumn:

  For `compareChromatograms()`: optional `character(1)` giving the name
  of a chromatogram variable (column in
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md))
  whose values should be used as row and column names of the returned
  array. The column must contain unique values. If `NULL` (the default),
  the array dimensions are unnamed. Only used when `y` is missing.

- threshold:

  For `peakBoundary()`: `numeric(1)` (default `0.1`). Fraction of the
  peak height above baseline used as a fallback cut-off when
  valley-based boundaries are not suitable. Must be `>= 0` and `< 1`.

- baselineThreshold:

  For `peakBoundary()`: `numeric(1)` (default `0.1`). Fraction of the
  peak height above the baseline. Valley positions returned by
  [`MsCoreUtils::valleys()`](https://rdrr.io/pkg/MsCoreUtils/man/valleys.html)
  are accepted only if the intensity at the valley is at or below
  `baseline + peak_height * baselineThreshold`. Must be `>= 0` and
  `< 1`.

- baselineQuantile:

  For `peakBoundary()`: `numeric(1)` (default `0.1`). Quantile of the
  intensity distribution used as the baseline estimate. Must be `>= 0`
  and `<= 1`.

## Value

Refer to the individual function description for information on the
return value.

## Filter Peaks Variables

Functions that filter a `Chromatograms`'s peaks data (i.e.,
`@peaksData`). These functions remove peaks data that do not meet the
specified conditions. If a chromatogram in a `Chromatograms` object is
filtered, only the corresponding peaks variable pairs (i.e., rows) in
the `peaksData` are removed, while the chromatogram itself remains in
the object.

The available functions to filter chromatographic peaks data include:

- `filterPeaksData()`: Filters numerical peaks data variables based on
  the specified numerical `ranges` parameter. This method returns the
  same input `Chromatograms` object, but the filtering step is added to
  the processing queue. The filtered data will be reflected when the
  user accesses `peaksData`. This function does *not* reduce the number
  of chromatograms in the object, but it removes the specified peaks
  data (e.g., "rtime" and "intensity" pairs) from the `peaksData`.

In the case of a read-only backend, (such as the
[ChromBackendMzR](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackendMzR.md)),
the replacement of the peaks data is not possible. The peaks data can be
filtered, but the filtered data will not be saved in the backend. This
means the original mzML files will not be affected by computations
performed on the
[Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md).

## Impute Peaks Variables

`imputePeaksData` will impute missing values in a `Chromatograms`'s
peaks data (i.e., `@peaksData`). This functions replace missing peaks
data values with specified imputation methods using various methods such
as linear interpolation, spline interpolation, Gaussian kernel
smoothing, or LOESS smoothing. This method modifies the peaks data in
place and returns the same `Chromatograms` object with imputed values.

## Peak Boundary Detection

`peakBoundary()` determines the retention time boundaries of the tallest
peak in each chromatogram. The function uses
[`MsCoreUtils::valleys()`](https://rdrr.io/pkg/MsCoreUtils/man/valleys.html)
to locate the valleys (local minima) flanking the apex. If the valley
intensities exceed a baseline-relative threshold (controlled by
`baselineThreshold`), it falls back to a threshold-based boundary search
using `threshold`. The baseline is estimated as the `baselineQuantile`
quantile of the chromatogram's intensity values. The result is a
`matrix` with one row per chromatogram and columns `left_boundary` and
`right_boundary` (retention times). Chromatograms that are empty, have
fewer than 3 data points, contain only `NA` or all-zero intensities
return `NA` for both boundaries.

## Compare Chromatograms

`compareChromatograms()` compares chromatograms in two steps:

1.  **Align** – `MAPFUN` (default `matchRtime()`) maps two chromatograms
    onto a common retention-time grid and returns `list(x, y)`, where
    `x` and `y` are numeric vectors of equal length containing the
    aligned intensities of the first and second chromatogram
    respectively.

2.  **Score** – `FUN` (default
    [`stats::cor()`](https://rdrr.io/r/stats/cor.html), Pearson
    correlation) computes a single similarity value from those aligned
    intensity vectors.

If `y` is missing, each chromatogram in `x` is compared against every
other chromatogram in `x`; otherwise, each in `x` is compared with each
in `y`.

The result is a 3-dimensional numeric array with dimensions `length(x)`
x `length(y)` x 2 (or symmetric `n x n x 2` for self-comparison). Layer
`[, , 1]` (named `"score"`) contains pairwise similarity scores; layer
`[, , 2]` (named `"n_peaks"`) contains the number of overlapping
retention-time points used to compute each score. Pairs with fewer
overlapping retention-time points than `minPeaks` (default 4) return
`NA` in the score layer; the actual overlap count is still recorded in
the `n_peaks` layer. The diagonal of a self-comparison is always 1
(score) and the number of data points in that chromatogram (count).

`matchRtime()` is the default `MAPFUN`. Given two chromatograms as
`data.frame`s with `rtime` and `intensity` columns, it aligns their RT
axes and returns a named list with elements `x` and `y`: equal-length
intensity vectors evaluated on a shared RT grid, ready for similarity
scoring.

The alignment works as follows: `matchRtime()` first identifies the RT
range where both chromatograms have measured points within `tolerance`
of each other (the *overlap*). Within that range, it builds a shared RT
grid from all of `x`'s RT points, adding any RT points from `y` that
have no close match in `x` (within `tolerance`). Both intensity vectors
are then linearly interpolated at grid positions they do not natively
cover, using
[`stats::approx()`](https://rdrr.io/r/stats/approxfun.html). If either
chromatogram has fewer than 2 data points, or the two chromatograms do
not overlap, empty vectors are returned.

The `tolerance` parameter (default `Inf`, meaning all RT points are
considered matching) controls the strictness of the matching. Lowering
it prevents comparing a measured peak against a long interpolated gap in
the other chromatogram. Pass `tolerance` via `...` in
`compareChromatograms()`.

When `y` is missing, the `labelsColumn` parameter assigns meaningful
row/column names to the output from a
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
column (e.g., `"mz"` or a user-defined feature identifier). The column
must contain unique values. To compare groups of chromatograms
separately, split the object with
[`split()`](https://rdrr.io/r/base/split.html) beforehand and apply
`compareChromatograms()` to each subset.

## See also

[Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
for a general description of the `Chromatograms` object, and
[chromData](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
for accessing,substituting and filtering chromatographic variables. For
more information on the queuing of processings and parallelization for
larger dataset processing see
[processingQueue](https://rformassspectrometry.github.io/Chromatograms/reference/processingQueue.md).

## Author

Philippine Louail

## Examples

``` r

# Create a Chromatograms object
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem2", "mem3")
)

pdata <- list(
    data.frame(
        rtime = c(12.4, 12.8, 13.2, 14.6),
        intensity = c(123.3, 153.6, 2354.3, 243.4)
    ),
    data.frame(
        rtime = c(45.1, 46.2),
        intensity = c(100, 80.1)
    ),
    data.frame(
        rtime = c(12.4, 12.8, 13.2, 14.6),
        intensity = c(123.3, 153.6, 2354.3, 243.4)
    )
)

be <- backendInitialize(new("ChromBackendMemory"),
    chromData = cdata,
    peaksData = pdata
)

chr <- Chromatograms(be)

# Access peaks data
peaksData(chr)
#> [[1]]
#>   rtime intensity
#> 1  12.4     123.3
#> 2  12.8     153.6
#> 3  13.2    2354.3
#> 4  14.6     243.4
#> 
#> [[2]]
#>   rtime intensity
#> 1  45.1     100.0
#> 2  46.2      80.1
#> 
#> [[3]]
#>   rtime intensity
#> 1  12.4     123.3
#> 2  12.8     153.6
#> 3  13.2    2354.3
#> 4  14.6     243.4
#> 

# Access specific peaks data variables
peaksData(chr, columns = "rtime")
#> [[1]]
#>   rtime
#> 1  12.4
#> 2  12.8
#> 3  13.2
#> 4  14.6
#> 
#> [[2]]
#>   rtime
#> 1  45.1
#> 2  46.2
#> 
#> [[3]]
#>   rtime
#> 1  12.4
#> 2  12.8
#> 3  13.2
#> 4  14.6
#> 
rtime(chr)
#> [[1]]
#> [1] 12.4 12.8 13.2 14.6
#> 
#> [[2]]
#> [1] 45.1 46.2
#> 
#> [[3]]
#> [1] 12.4 12.8 13.2 14.6
#> 

# Replace peaks data
rtime(chr)[[1]] <- c(1, 2, 3, 4)

# Filter peaks data
filterPeaksData(chr, variables = "rtime", ranges = c(12.5, 13.5))
#> Chromatographic data (Chromatograms) with 3 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       1 134.4
#> ... 3 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> Lazy evaluation queue: 1 processing step(s)
#> Processing:
#>  Filter: remove peaks based on the variables: rtimethe ranges: 12.5, 13.5and the match condition: any [Wed Aug  5 11:48:06 2026]
#>  Filter: remove peaks based on the variables: rtimethe ranges: 12.5, 13.5and the match condition: all [Wed Aug  5 11:48:06 2026] 

# Pairwise similarity: returns a 3D array [i, j, layer]
res <- compareChromatograms(chr)
res[, , "score"]   ## similarity scores
#>      [,1] [,2] [,3]
#> [1,]    1   NA   NA
#> [2,]   NA    1   NA
#> [3,]   NA   NA    1
res[, , "n_peaks"] ## number of overlapping RT points
#>      [,1] [,2] [,3]
#> [1,]    4    0    0
#> [2,]    0    2    0
#> [3,]    0    0    4

## Use Spearman correlation (passed to cor() via ...)
compareChromatograms(chr, method = "spearman")[, , "score"]
#>      [,1] [,2] [,3]
#> [1,]    1   NA   NA
#> [2,]   NA    1   NA
#> [3,]   NA   NA    1

# Use a chromData column as row/column labels
compareChromatograms(chr, labelsColumn = "mz")[, , "score"]
#>       112.2 123.3 134.4
#> 112.2     1    NA    NA
#> 123.3    NA     1    NA
#> 134.4    NA    NA     1

# Compare two Chromatograms objects
compareChromatograms(chr[1:2], chr[3])
#> , , score
#> 
#>      [,1]
#> [1,]   NA
#> [2,]   NA
#> 
#> , , n_peaks
#> 
#>      [,1]
#> [1,]    0
#> [2,]    0
#> 
```
