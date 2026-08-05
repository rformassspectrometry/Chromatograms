# Chromatographic MS Data Backends

`ChromBackend` is a virtual class that defines what different backends
need to provide to be used by the `Chromatograms` package and classes.

The backend should provide access to the chromatographic data which
mainly consists of (paired) intensity and retention time values.
Additional chromatographic metadata such as MS level and precursor and
product m/z should also be provided.

Through their implementation different backends can be either optimized
for minimal memory requirements or performance. Each backend needs to
implement data access methods listed in section *Backend functions:*
below.

And example implementation and more details and descriptions are
provided in the *Creating new `ChromBackend` classes for Chromatograms*
vignette.

Currently available backends are:

- `ChromBackendMemory`: This backend stores chromatographic data
  directly in memory, making it ideal for small datasets or testing. It
  can be initialized with a `data.frame` of chromatographic data via the
  `chromData` parameter and a `list` of `data.frame` entries for peaks
  data using the `peaksData` parameter. These data can be accessed with
  the
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  and
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  functions.

- `ChromBackendMzR`: The `ChromBackendMzR` inherits all slots and
  methods from the base `ChromBackendMemory` backend, providing
  additional functionality for reading chromatographic data from mzML
  files.

- `ChromBackendSpectra`: The `ChromBackendSpectra` inherits all slots
  and methods from the base `ChromBackendMemory` backend, providing
  additional functionality for reading chromatographic data from
  `Spectra` objects.

Filter the peak data based on the provided ranges for the given
variables.

## Usage

``` r
coreChromVariables()

corePeaksVariables()

# S4 method for class 'ChromBackend'
x$name

# S4 method for class 'ChromBackend'
x$name <- value

# S4 method for class 'ChromBackend'
backendMerge(object, ...)

# S4 method for class 'ChromBackend'
chromData(object, columns = chromVariables(object), drop = FALSE)

# S4 method for class 'ChromBackend'
chromData(object) <- value

# S4 method for class 'ChromBackend'
chromExtract(object, peak.table, by)

# S4 method for class 'ChromBackend'
factorize(object, ...)

# S4 method for class 'ChromBackend'
peaksData(object, columns = c("rtime", "intensity"), drop = FALSE, ...)

# S4 method for class 'ChromBackend'
peaksData(object) <- value

# S4 method for class 'ChromBackend'
x[i, j, ..., drop = FALSE]

# S4 method for class 'ChromBackend'
x[[i, j, ...]]

# S4 method for class 'ChromBackend'
x[[i, j, ...]] <- value

# S4 method for class 'ChromBackend'
backendBpparam(object, BPPARAM = bpparam())

# S4 method for class 'ChromBackend'
backendInitialize(object, ...)

# S4 method for class 'ChromBackend'
backendParallelFactor(object, ...)

# S4 method for class 'list'
backendMerge(object, ...)

# S4 method for class 'ChromBackend'
chromIndex(object)

# S4 method for class 'ChromBackend'
chromIndex(object) <- value

# S4 method for class 'ChromBackend'
chromVariables(object)

# S4 method for class 'ChromBackend'
collisionEnergy(object)

# S4 method for class 'ChromBackend'
collisionEnergy(object) <- value

# S4 method for class 'ChromBackend'
dataOrigin(object)

# S4 method for class 'ChromBackend'
dataOrigin(object) <- value

# S4 method for class 'ChromBackend,ANY'
extractByIndex(object, i)

# S4 method for class 'ChromBackend,missing'
extractByIndex(object, i)

# S4 method for class 'ChromBackend'
intensity(object)

# S4 method for class 'ChromBackend'
intensity(object) <- value

# S4 method for class 'ChromBackend'
isEmpty(x)

# S4 method for class 'ChromBackend'
isReadOnly(object)

# S4 method for class 'ChromBackend'
length(x)

# S4 method for class 'ChromBackend'
lengths(x)

# S4 method for class 'ChromBackend'
msLevel(object)

# S4 method for class 'ChromBackend'
msLevel(object) <- value

# S4 method for class 'ChromBackend'
mz(object)

# S4 method for class 'ChromBackend'
mz(object) <- value

# S4 method for class 'ChromBackend'
mzMax(object)

# S4 method for class 'ChromBackend'
mzMax(object) <- value

# S4 method for class 'ChromBackend'
mzMin(object)

# S4 method for class 'ChromBackend'
mzMin(object) <- value

# S4 method for class 'ChromBackend'
peaksVariables(object)

# S4 method for class 'ChromBackend'
precursorMz(object)

# S4 method for class 'ChromBackend'
precursorMz(object) <- value

# S4 method for class 'ChromBackend'
precursorMzMax(object)

# S4 method for class 'ChromBackend'
precursorMzMax(object) <- value

# S4 method for class 'ChromBackend'
precursorMzMin(object)

# S4 method for class 'ChromBackend'
precursorMzMin(object) <- value

# S4 method for class 'ChromBackend'
productMz(object)

# S4 method for class 'ChromBackend'
productMz(object) <- value

# S4 method for class 'ChromBackend'
productMzMax(object)

# S4 method for class 'ChromBackend'
productMzMax(object) <- value

# S4 method for class 'ChromBackend'
productMzMin(object)

# S4 method for class 'ChromBackend'
productMzMin(object) <- value

# S4 method for class 'ChromBackend'
reset(object)

# S4 method for class 'ChromBackend'
rtime(object)

# S4 method for class 'ChromBackend'
rtime(object) <- value

# S4 method for class 'ChromBackend,ANY'
split(x, f, drop = FALSE, ...)

# S4 method for class 'ChromBackend'
filterChromData(
  object,
  variables = character(),
  ranges = numeric(),
  match = c("any", "all"),
  keep = TRUE
)

# S4 method for class 'ChromBackend'
filterEmptyChromatograms(object, ...)

# S4 method for class 'ChromBackend'
filterPeaksData(
  object,
  variables = character(),
  ranges = numeric(),
  match = c("any", "all"),
  keep = TRUE
)

# S4 method for class 'ChromBackend'
supportsSetBackend(object, ...)

# S4 method for class 'ChromBackend'
imputePeaksData(
  object,
  method = c("linear", "spline", "gaussian", "loess"),
  span = 0.3,
  sd = 1,
  window = 2,
  extrapolate = FALSE,
  ...
)
```

## Arguments

- x:

  Object extending `ChromBackend`.

- name:

  For `$` and `$<-`: the name of the chromatogram variable to return or
  set.

- value:

  Replacement value for `<-` methods. See individual method description
  or expected data type.

- object:

  Object extending `ChromBackend`.

- ...:

  Additional arguments.

- columns:

  For
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  accessor: optional `character` with column names (chromatogram
  variables) that should be included in the returned `data.frame`. By
  default, all columns are returned.

- drop:

  For
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  and
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  `logical(1)` default to `FALSE`. If `TRUE`, and one column is
  requested by the user, the method should return a vector (or list of
  vector for
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md))
  of the single column requested.

- peak.table:

  For
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  A `data frame` containing the following minimum columns:

  - rtMin: Minimum retention time for each peak. Cannot be NA.

  - rtMax: Maximum retention time for each peak. Cannot be NA.

  - mzMin: Minimum m/z value for each peak.

  - mzMax: Maximum m/z value for each peak. Additionally, the
    `peak.table` must include columns that uniquely identify
    chromatograms in the `object`. Common choices are "msLevel" and/or
    "dataOrigin". These columns must also be present in the `chromData`
    of the `object`. Any extra columns in `peak.table` will be added to
    the `chromData` of the newly created object.

- by:

  for
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)` `A
  `character` vector specifying one or more column names that are
  present in both `peak.table` and `chromData(object)`. These columns
  uniquely identify chromatograms. The combination of these columns must
  be unique in `chromData(object)`. Can be of length 1 or greater.

- i:

  For `[`: `integer`, `logical` or `character` to subset the object.

- j:

  For `[` and `[[`: ignored.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- f:

  `factor` defining the grouping to split `x`. See
  [`split()`](https://rdrr.io/r/base/split.html).

- variables:

  For
  [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  `character` vector with the names of the chromatogram variables to
  filter for. The list of available chromatogram variables can be
  obtained with
  [`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md).

- ranges:

  For
  [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  : a `numeric` vector of paired values (upper and lower boundary) that
  define the ranges to filter the `object`. These paired values need to
  be in the same order as the `variables` parameter (see below).

- match:

  For
  [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  : `character(1) ` defining whether the condition has to match for all
  provided `ranges` (`match = "all"`; the default), or for any of them
  (`match = "any"`) for chromatogram data to be retained.

- keep:

  For
  [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  `logical(1)` defining whether to keep (`keep = TRUE`) or remove
  (`keep = FALSE`) the chromatogram data that match the condition.

- method:

  For
  [`imputePeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  `character(1)`: Imputation method ("linear", "spline", "gaussian",
  "loess")

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

## Value

Refer to the individual function description for information on the
return value.

## Core chromatogram variables

The *core* chromatogram variables are variables (metadata) that
can/should be provided by a backend. For each of these variables a value
needs to be returned, if none is defined, a missing value (of the
correct data type) should be returned. The names of the chromatogram
variables in your current chromatogram object are returned with the
[`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
function.

For each core chromatogram variable a dedicated access method exists. In
contrast to the peaks data described below, a single value should be
returned for each chromatogram.

The `coreChromVariables()` function returns the core chromatogram
variables along with their expected (defined) data type.

The core chromatogram variables (in alphabetical order) are:

- `chromIndex`: an `integer` with the index of the chromatogram in the
  original source file (e.g. *mzML* file). In backedn with no original
  source file, this variable should be set to `NA_integer_`.

- `collisionEnergy`: for SRM data, `numeric` with the collision energy
  of the precursor.

- `dataOrigin`: optional `character` with the origin of a chromatogram.

- `msLevel`: `integer` defining the MS level of the data.

- `mz`: optional `numeric` with the (target) m/z value for the
  chromatographic data.

- `mzMin`: optional `numeric` with the lower m/z value of the m/z range
  in case the data (e.g. an extracted ion chromatogram EIC) was
  extracted from a `Spectra` object.

- `mzMax`: optional `numeric` with the upper m/z value of the m/z range.

- `precursorMz`: for SRM data, `numeric` with the target m/z of the
  precursor (parent).

- `precursorMzMin`: for SRM data, optional `numeric` with the lower m/z
  of the precursor's isolation window.

- `precursorMzMax`: for SRM data, optional `numeric` with the upper m/z
  of the precursor's isolation window.

- `productMz` for SRM data, `numeric` with the target m/z of the product
  ion.

- `productMzMin`: for SRM data, optional `numeric` with the lower m/z of
  the product's isolation window.

- `productMzMax`: for SRM data, optional `numeric` with the upper m/z of
  the product's isolation window.

## Core Peaks variables

Similar to the *core* chromatogram variables, *core* peaks variables
represent metadata that should be provided by a backend. Each of these
variables should return a value, and if undefined, a missing value (with
the appropriate data type) is returned. The number of values for a peaks
variable in a single chromatogram can vary, from none to multiple, and
may differ between chromatograms.

The names of peaks variables in the current chromatogram object can be
obtained with the
[`peaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
function.

Each core peaks variable has a dedicated accessor method.

The `corePeaksVariables()` function returns the core peaks variables
along with their expected (defined) data type.

The core peaks variables, listed in the required order for `peaksData`,
are:

- `rtime`: A `numeric` vector containing retention time values.

- `intensity`: A `numeric` vector containing intensity values.

They should be provided for each chromatogram in the backend, **in this
order**, No NAs are allowed for the `rtime` values. These
characteristics will be checked with the
[`validPeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
function.

## Mandatory methods

New backend classes **must** extend the base `ChromBackend` class and
implement the following mandatory methods:

- `backendInitialize()`: initialises the backend. This method is
  supposed to be called right after creating an instance of the backend
  class and should prepare the backend. Parameters can be defined freely
  for each backend, depending on what is needed to initialize the
  backend. This method has to ensure to set the chromatogram variable
  `dataOrigin` correctly.

- `backendBpparam()`: returns the parallel processing setup supported by
  the backend class. This function can be used by any higher level
  function to evaluate whether the provided parallel processing setup
  (or the default one returned by `bpparam()`) is supported by the
  backend. Backends not supporting parallel processing (e.g. because
  they contain a connection to a database that can not be shared across
  processes) should extend this method to return only `SerialParam()`
  and hence disable parallel processing for (most) methods and
  functions. See also `backendParallelFactor()` for a function to
  provide a preferred splitting of the backend for parallel processing.

- `backendParallelFactor()`: returns a `factor` defining an optimal
  (preferred) way how the backend can be split for parallel processing
  used for all peak data accessor or data manipulation functions. The
  default implementation returns a factor of length 0
  ([`factor()`](https://rdrr.io/r/base/factor.html)) providing thus no
  default splitting. `backendParallelFactor()` for `ChromBackendMzR` on
  the other hand returns `factor(dataOrigin(object))` hence suggesting
  to split the object by data file.

- [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),
  `chromData<-`: gets or sets general chromatogram metadata
  (annotation).
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `data.frame`, `chromData<-` expects a `data.frame` with the
  same number of rows as there are chromatograms in `object`. Read-only
  backends might not need to implement the replacement method
  `chromData<-` (unless some internal caching mechanism could be used).
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  should be implemented with the parameter `drop` set to `FALSE` as
  default. With `drop = FALSE` the method should return a `data.frame`
  even if one column is requested. If `drop = TRUE` is specified, the
  output will be a vector of the single column requested. New backends
  should be implemented such as if empty, the method returns a
  `data.frame` with 0 rows and the columns defined by
  [`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md).
  By default, the function *should* return at minimum the
  coreChromVariables, even if NAs.

- [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md):
  return A new `Chrombackend` object containing separated
  chromatographic area as individual chromatograms. The chromatographic
  areas are defined by the `peak.table` parameter. The new object will
  contain chromatograms that match the conditions defined in
  `peak.table`. If no chromatograms match the conditions, an empty
  `ChromBackend` object should be returned.

- `extractByIndex()`: function to subset a backend to selected elements
  defined by the provided index. Similar to `[`, this method should
  allow extracting (or to subset) the data in any order. In contrast to
  `[`, however, `i` is expected to be an `integer` (while `[` should
  also support `logical` and eventually `character`). While being
  apparently redundant to `[`, this methods avoids package namespace
  errors/problems that can result in implementations of `[` being not
  found by R (which can happen sometimes in parallel processing using
  the
  [`BiocParallel::SnowParam()`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html)).
  This method is used internally to extract/subset its backend.
  Implementation of this method is mandatory.

- [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  returns a `list` of `data.frame` with the data (e.g. retention time -
  intensity pairs) from each chromatogram. The length of the `list` is
  equal to the number of chromatograms in `object`. For an empty
  chromatogram a `data.frame` with 0 rows and two columns (named
  `"rtime"` and `"intensity"`) has to be returned. The optional
  parameter `columns`, if supported by the backend allows to define
  which peak variables should be returned in each array. As default
  (minimum) columns `"rtime"` and `"intensity"` have to be provided.
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  should be implemented with the parameter `drop` set to `FALSE` as
  default. With `drop = FALSE` the method should return a `data.frame`
  even if only one column is requested. If `drop = TRUE` is specified,
  the output will be a vector of the single column requested.

- `peaksData<-` replaces the peak data (retention time and intensity
  values) of the backend. This method expects a `list` of
  two-dimensional arrays (`data.frame`) with columns representing the
  peak variables. All existing peaks data are expected to be replaced
  with these new values. The length of the `list` has to match the
  number of chromatogram of `object`. Note that only writeable backends
  need to support this method.

- `[`: subset the backend. Only subsetting by element (*row*/`i`) is
  allowed. This method should be implemented as to support empty
  integer.

- `$`, `$<-`: access or set/add a single chromatogram variable (column)
  in the backend.

- `backendMerge()`: merges (combines) `ChromBackend` objects into a
  single instance. All objects to be merged have to be of the same type.

## Optional methods with default implementations

Additional methods that might be implemented, but for which default
implementations are already present are:

- `[[`

- `backendParallelFactor()`: returns a `factor` defining an optimal
  (preferred) way how the backend can be split for parallel processing
  used for all *peak* data accessor or data manipulation functions. The
  default implementation returns a factor of length 0
  ([`factor()`](https://rdrr.io/r/base/factor.html)) providing thus no
  default splitting.

- [`chromIndex()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  returns an `integer` vector with the index of the chromatograms in the
  original source file.

- [`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  returns a `character` vector with the available chromatogram variables
  (columns, fields or attributes) available in `object`. Variables
  listed by this function are expected to be returned (if requested) by
  the
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  function.

- [`collisionEnergy()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),
  `collisionEnergy<-`: gets or sets the collision energy for the
  precursor (for SRM data).
  [`collisionEnergy()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  `object`.

- [`dataOrigin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),
  `dataOrigin<-`: gets or sets the *data origin* variable.
  [`dataOrigin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `character` of length equal to the number of chromatograms,
  `dataOrigin<-` expects a `character` of length equal `length(object)`.

- [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  filters any numerical chromatographic data variables based on the
  provided numerical `ranges`. The method should return a `ChromBackend`
  object with the chromatograms that match the condition. This function
  will results in an object with less chromatogram than the original.

- [`filterEmptyChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md):
  removes empty chromatograms (i.e. chromatograms without peaks).
  Implementation of this method is optional since a default
  implementation for `ChromBackend` is available.

- `intensity()`: gets the intensity values from the chromatograms.
  Returns a `list` of `numeric` vectors (intensity values for each
  chromatogram). The length of the list is equal to the number of
  chromatograms in `object`.

- `intensity<-`: replaces the intensity values. `value` has to be a
  `list` of length equal to the number of chromatograms and the number
  of values within each list element identical to the number of data
  pairs in each chromatogram. Note that just writeable backends need to
  support this method.

- [`imputePeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  Imputes missing intensity values in the chromatographic peaks data
  using various methods such as linear interpolation, spline
  interpolation, Gaussian kernel smoothing, or LOESS smoothing. This
  method modifies the peaks data in place and returns the same
  `ChromBackend` object with imputed values.

- `isReadOnly()`: returns a `logical(1)` whether the backend is *read
  only* or does allow also to write/update data. Defaults to FALSE.

- `isEmpty()`: returns a `logical` of length equal to the number of
  chromatograms with `TRUE` for chromatograms without any data pairs.

- [`length()`](https://rdrr.io/r/base/length.html): returns the number
  of chromatograms in the object.

- [`lengths()`](https://rdrr.io/r/base/lengths.html): returns the number
  of data pairs (retention time and intensity values) per chromatogram.

- [`msLevel()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  gets the chromatogram's MS level. Returns an `integer` vector (of
  length equal to the number of chromatograms) with the MS level for
  each chromatogram (or `NA_integer_` if not available).

- [`mz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),`mz<-`:
  gets or sets the m/z value of the chromatograms.
  [`mz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  ` object`, `mz<-` expects a `numeric` of length `length(object)`.

- [`mzMax()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),`mzMax<-`:
  gets or sets the upper m/z of the mass-to-charge range from which a
  chromatogram contains signal (e.g. if the chromatogram was extracted
  from MS data in spectra format and a m/z range was provided).
  [`mzMax()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  `object`, `mzMax<-` expects a `numeric` of length equal to the number
  of chromatograms in `object`.

- [`mzMin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),`mzMin<-`:
  gets or sets the lower m/z of the mass-to-charge range from which a
  chromatogram contains signal (e.g. if the chromatogram was extracted
  from MS data in spectra format and a m/z range was provided).
  [`mzMin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  `object`, `mzMin<-` expects a `numeric` of length equal to the number
  of chromatograms in `object`.

- [`peaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  lists the available data variables for the chromatograms. Default peak
  variables are `"rtime"` and `"intensity"` (which all backends need to
  support and provide), but some backends might provide additional
  variables. Variables listed by this function are expected to be
  returned (if requested) by the
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  function.

- [`precursorMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),`precursorMz<-`:
  gets or sets the (target) m/z of the precursor (for SRM data).
  [`precursorMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  `object`. `precursorMz<-` expects a `numeric` of length equal to the
  number of chromatograms.

- [`precursorMzMin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),[`precursorMzMax()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),[`productMzMin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),
  [`productMzMax()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md):
  gets the lower and upper margin for the precursor or product isolation
  windows. These functions might return the value of
  [`productMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  if the respective minimal or maximal m/z values are not defined in
  `object`.

- [`productMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md),`productMz<-`:
  gets or sets the (target) m/z of the product (for SRM data).
  [`productMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  returns a `numeric` of length equal to the number of chromatograms in
  `object`. `productMz<-` expects a `numeric` of length equal to the
  number of chromatograms.

- `rtime()`: gets the retention times from the chromatograms. returns a
  `list` of `numeric` vectors (retention times for each chromatogram).
  The length of the returned list is equal to the number of
  chromatograms in `object`.

- `rtime<-`: replaces the retention times. `value` has to be a `list` of
  length equal to the number of chromatograms and the number of values
  within each list element identical to the number of data pairs in each
  chromatogram. Note that just writeable backends support this method.

- [`split()`](https://rdrr.io/r/base/split.html): splits the backend
  into a `list` of backends (depending on parameter `f`). The default
  method for `ChromBackend` uses
  [`split.default()`](https://rdrr.io/r/base/split.html), thus backends
  extending `ChromBackend` don't necessarily need to implement this
  method.

- `supportsSetBackend()`: whether a `ChromBackend` supports the
  `Chromatograms` `setBackend()` function. The default function will
  take the
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  and
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  of the user's backend and pass it to the new backend. If the backend
  does not support this function, it should return `FALSE`. Therefore
  both backend in question should have a adequate
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  and
  [`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  method as well as their respective replacement method.

## Implementation notes

Backends extending `ChromBackend` **must** implement all of its methods
(listed above). A guide to create new backend classes is provided as a
dedicated vignette. Additional information and an example for a backend
implementation is provided in the respective vignette.

## Author

Johannes Rainer, Philippine Louail

## Examples

``` r

## Create a simple backend implementation
ChromBackendDummy <- setClass("ChromBackendDummy",
    contains = "ChromBackend"
)

## We will show examples on a `ChromBackendMemory` backend.
be <- ChromBackendMemory()

## The `backendInitialize()` method initializes the backend filling it with
## data. This method can take any parameters needed for the backend to
## get loaded with the data.
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
be <- backendInitialize(be, chromData = cdata, peaksData = pdata)

be
#> ChromBackendMemory with 3 chromatograms
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       1 134.4
#> ... 3 more  chromatogram variables/columns
#> ... 2 peaksData variables

## Data can be accessed with the accessor methods
msLevel(be)
#> [1] 1 1 1

rtime(be)
#> [[1]]
#> [1] 12.4 12.8 13.2 14.6
#> 
#> [[2]]
#> [1] 45.1 46.2
#> 
#> [[3]]
#> [1] 12.4 12.8 13.2 14.6
#> 

## Even if no data was provided for all chromatogram variables, its accessor
## methods are supposed to return a value.
precursorMz(be)
#> [1] NA NA NA

## The `peaksData()` method is supposed to return data/frames of rtime and
## intensity pairs as a `list`.
peaksData(be)
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

## Use columns to extract specific peaks variables. Below we extract rtime
## and intensity values, but in reversed order to the default.
peaksData(be, columns = c("intensity", "rtime"))
#> [[1]]
#>   intensity rtime
#> 1     123.3  12.4
#> 2     153.6  12.8
#> 3    2354.3  13.2
#> 4     243.4  14.6
#> 
#> [[2]]
#>   intensity rtime
#> 1     100.0  45.1
#> 2      80.1  46.2
#> 
#> [[3]]
#>   intensity rtime
#> 1     123.3  12.4
#> 2     153.6  12.8
#> 3    2354.3  13.2
#> 4     243.4  14.6
#> 

## List available chromatographic variables
chromVariables(be)
#>  [1] "msLevel"         "mz"              "dataOrigin"      "chromIndex"     
#>  [5] "collisionEnergy" "mzMin"           "mzMax"           "precursorMz"    
#>  [9] "precursorMzMin"  "precursorMzMax"  "productMz"       "productMzMin"   
#> [13] "productMzMax"   

## List available peak variables
peaksVariables(be)
#> [1] "rtime"     "intensity"

## Extract multiple chromatographic variables
chromData(be, c("dataOrigin", "mz", "msLevel"))
#>   dataOrigin    mz msLevel
#> 1       mem1 112.2       1
#> 2       mem2 123.3       1
#> 3       mem3 134.4       1

## Single variables can also be accessed and replaced
mz(be)
#> [1] 112.2 123.3 134.4
mz(be) <- c(123.4, 134.5, 145.6)

be$msLevel
#> [1] 1 1 1
be$msLevel <- c(2L, 2L, 2L)

be[["rtime"]]
#> [[1]]
#> [1] 12.4 12.8 13.2 14.6
#> 
#> [[2]]
#> [1] 45.1 46.2
#> 
#> [[3]]
#> [1] 12.4 12.8 13.2 14.6
#> 
be[["rtime"]] <- list(
    c(12.4, 12.8, 13.2, 14.6),
    c(45.1, 46.2),
    c(12.4, 12.8, 13.2, 14.6)
)
```
