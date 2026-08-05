# The Chromatograms class to manage and access chromatographic data

The `Chromatograms` class encapsules chromatographic data and related
metadata. The chromatographic data is represented by a *backend*
extending the virtual
[ChromBackend](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
class which provides the raw data to the `Chromatograms` object.
Different backends and their properties are described in the
[ChromBackend](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
class documentation.

**Available Backends**: The package provides several backends:

- `ChromBackendMemory`: Stores data in memory (default, ideal for small
  datasets).

- `ChromBackendMzR`: Reads peaks data from raw MS files on demand.

- `ChromBackendSpectra`: Generates chromatographic data from a `Spectra`
  object. This backend supports both in-memory and file-backed `Spectra`
  objects, using an internal `spectraSortIndex` to avoid physically
  reordering the spectra.

## Usage

``` r
# S4 method for class 'ChromBackendOrMissing'
Chromatograms(object = ChromBackendMemory(), processingQueue = list(), ...)

# S4 method for class 'Spectra'
Chromatograms(
  object,
  summarize.method = c("sum", "max"),
  chromData = data.frame(),
  factorize.by = c("msLevel", "dataOrigin"),
  spectraVariables = character(),
  ...
)

# S4 method for class 'Chromatograms,ChromBackend'
setBackend(
  object,
  backend,
  f = processingChunkFactor(object),
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'Chromatograms'
x$name

# S4 method for class 'Chromatograms'
x$name <- value

# S4 method for class 'Chromatograms'
x[i, j, ..., drop = FALSE]

# S4 method for class 'Chromatograms'
x[[i, j, ...]]

# S4 method for class 'Chromatograms'
x[[i, j, ...]] <- value

# S4 method for class 'Chromatograms'
factorize(object, factorize.by = c("msLevel", "dataOrigin"), ...)

# S4 method for class 'Chromatograms'
chromExtract(object, peak.table, by, ...)

# S4 method for class 'Chromatograms'
filterEmptyChromatograms(object, ...)
```

## Arguments

- object:

  A Chromatograms object.

- processingQueue:

  [list](https://rdrr.io/r/base/list.html) a list of processing steps
  (i.e. functions) to be applied to the chromatographic data. The
  processing steps are applied in the order they are listed in the
  `processingQueue`.

- ...:

  Additional arguments.

- summarize.method:

  For `Chromatograms` created with a `Spectra` object: A `character`
  vector with the name of the function to be used to summaries the
  spectra data intensity. The available methods are "sum" and "max". The
  default is "sum".

- chromData:

  For `Chromatograms()` build from a `Spectra` object backend, a
  `data.frame` with the chromatographic data. If not provided (or if
  empty), a default `data.frame` with the core chromatographic variables
  will be created.

- factorize.by:

  A `character` vector with the names of the variables in the `Spectra`
  object and the `chromData` slot that should be used to factorize the
  `Spectra` object data to generate the chromatographic data.

- spectraVariables:

  A `character` vector specifying which variables from the `Spectra`
  object should be added to the chromData. These will be mapped using
  the `chromSpectraIndex` variable.

- backend:

  [ChromBackend](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
  object providing the raw data for the `Chromatograms` object.

- f:

  `factor` defining the grouping to split the `Chromatograms` object.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- x:

  A Chromatograms object.

- name:

  A `character` string specifying the name of the variable to access.

- value:

  The value to replace the variable with.

- i:

  For `[`: `integer`, `logical` or `character` to subset the object.

- j:

  For `[` and `[[`: ignored.

- drop:

  For `[`: `logical(1)` default to `FALSE`.

- peak.table:

  For `chromExtract()` A `data frame` containing the following minimum
  columns:

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

  A `character` vector naming one or more columns that uniquely identify
  chromatograms in both `peak.table` and `chromData(object)`. The
  combination of these columns must be unique within
  `chromData(object)`. Typically includes `"dataOrigin"`, `"msLevel"`,
  or both.

## Value

Refer to the individual function description for information on the
return value.

## Note

This needs to be discussed, if we want for example to be able to set a a
backend to `ChromBackendMzR` we need to implement backendInitialize()
better. = Support peaksData and chromData as arguments AND have a way to
write .mzml files (which we do not have for chromatographic data).

## Creation of objects

`Chromatograms` objects can be created using the `Chromatograms()`
construction function. Either by providing a `ChromBackend` object or by
providing a `Spectra` object. The `Spectra` object will be used to
generate a `Chromatograms` object with a backend of class
[`ChromBackendSpectra`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackendSpectra.md).

## Data stored in a `Chromatograms` object

The `Chromatograms` object is a container for chromatographic data,
which includes peaks data (*retention time* and related intensity
values, also referred to as *peaks data variables* in the context of
`Chromatograms`) and metadata of individual chromatogram (so called
*chromatograms variables*). While a core set of chromatograms variables
(the `coreChromatogramsVariables()`) and peaks data variables (the
[`corePeaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md))
are guaranteed to be provided by a `Chromatograms`, it is possible to
add arbitrary variables to a `Chromatograms` object.

The `Chromatograms` object is designed to contain chromatographic data
of a (large) set of chromatograms. The data is organized *linearly* and
can be thought of a list of chromatograms, i.e. each element in the
`Chromatograms` is one chromatogram.

The *chromatograms variables* information in the `Chromatograms` object
can be accessed using the
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
function. Specific chromatograms variables can be accessed by either
precising the `"columns"` parameter in
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
or using `$`. `@chromData` can be accessed, replaced but also
filtered/subsetted. Refer to the
[chromData](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
documentation for more details.

The *peaks data variables* information in the `Chromatograms` object can
be accessed using the
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
function. Specific peaks variables can be accessed by either precising
the `"columns"` parameter in
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
or using `$`. `@peaksData` can be accessed, replaced but also
filtered/subsetted. Refer to the
[peaksData](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
documentation for more details.

## Processing of `Chromatograms` objects

Functions that process the chromatograms data in some ways can be
applied to the object either directly or by using the `processingQueue`
mechanism. The `@processingQueue` is a list of processing steps that are
stored within the object and only applied when needed. This was created
so that the data can be processed in a single step and is very useful
for larger datasets. This is even more true as this processing queue
will call function that can be applied on the data in a chunk-wise
manner. This allows for parallel processing of the data and reduces the
memory demand. To read more about the `processingQueue`, and how to
parallelize your processes, see the
[processingQueue](https://rformassspectrometry.github.io/Chromatograms/reference/processingQueue.md)
documentation.

## Subsetting and accessing data

The `Chromatograms` class supports subsetting by chromatogram (i.e.
rows) using the `[` operator. The `[` operator does not support
subsetting by columns. Specific chromatograms or peaks variables can be
accessed using the `[[` operator or the `$` operator. The `[[` operator
can also be used to replace specific chromatograms or peaks variables.

## Changing the backend

The
[`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function can be used to change the backend of a `Chromatograms` object.
This can be useful to switch to a backend that better suits the needs of
the user, for example switching to a memory-based backend for smaller
datasets or to a file-based backend for larger datasets. The
[`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function supports parallelization of the backend conversion using the
`BPPARAM` parameter. Note that any queued processing steps are applied
during the backend switch (since
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
is called to transfer the data) and the processing queue is emptied
afterwards.

## Filtering chromatograms

- `filterEmptyChromatograms()`: removes empty chromatograms (i.e.
  chromatograms without peaks) from the object. Returns the filtered
  `Chromatograms` object (with chromatograms in their original order).

## Extracting chromatograms based on a peak table

The `chromExtract()` function allows users to extract specific regions
of interest from a `Chromatograms` object based on a user-provided peak
table. Each row in the `peak.table` defines a region to extract, using
minimum and maximum retention time (and m/z in the case of
`chromBackendSpectra`) boundaries, and identifiers that uniquely match
chromatograms in the object.

The resulting **new** `Chromatograms` object contains only chromatograms
overlapping the specified regions, with updated metadata reflecting the
extracted boundaries.

This function is most commonly used to subset chromatographic data
around detected peaks or predefined time/mass ranges, for example to
reprocess, visualize, or quantify extracted chromatograms corresponding
to known features. It's important to notes that filtering by m/z is only
supported when using a `ChromBackendSpectra` backend. if the `mzMin` and
`mzMax` columns are provided when using other backends, they will be
ignored.

## See also

[chromData](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
for a general description of the chromatographic metadata available in
the object, as well as how to access, replace and subset them.
[peaksData](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
for a general description of the chromatographic peaks data available in
the object, as well as how to access, replace and subset them.
[processingQueue](https://rformassspectrometry.github.io/Chromatograms/reference/processingQueue.md)
for more information on the queuing of processings and parallelization
for larger dataset.

## Examples

``` r

## Create a Chromatograms object with ChromBackendMemory
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem2", "mem3")
)
pdata <- list(
    data.frame(rtime = c(2.1, 2.5, 3.0, 3.4, 3.9),
               intensity = c(100, 250, 400, 300, 150)),
    data.frame(rtime = c(3.5, 4.0, 4.5),
               intensity = c(80, 120, 90)),
    data.frame(rtime = c(5.1, 5.8, 6.3, 6.9, 7.5),
               intensity = c(80, 500, 1200, 600, 120))
)
chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
chr
#> Chromatographic data (Chromatograms) with 3 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       1 134.4
#> ... 3 more  chromatogram variables/columns
#> ... 2 peaksData variables

## Create a Chromatograms object from a Spectra object
library(MsBackendMetaboLights)
library(Spectra)

be <- backendInitialize(MsBackendMetaboLights(),
    mtblsId = "MTBLS39",
    filePattern = c("63B.cdf")
)
#> Used data files from the assay's column "Raw Spectral Data File" since none were available in column "Derived Spectral Data File".
s <- Spectra(be)
s <- setBackend(s, MsBackendMemory())
chr <- Chromatograms(s)

## Subset
chr[1:2]
#> Chromatographic data (Chromatograms) with 2 chromatograms in a ChromBackendSpectra backend:
#>   chromIndex msLevel  mz
#> 1         NA       1 Inf
#> 2         NA       1 Inf
#> ... 7 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> 
#> The Spectra object contains 1101 spectra

## Access a specific variable
chr[["msLevel"]]
#> [1] 1 1 1
chr$msLevel
#> [1] 1 1 1

## Replace data of a specific variable
chr$msLevel <- c(2L, 2L, 2L)

## Re-factorize the data
chr <- factorize(chr)

## Change the backend to memory
chr <- setBackend(chr, ChromBackendMemory())
```
