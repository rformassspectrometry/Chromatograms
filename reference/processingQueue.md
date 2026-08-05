# Efficiently processing `Chromatograms` objects.

The `processingQueue` of a `Chromatograms` object is a list of
processing steps (i.e., functions) that are stored within the object and
applied only when needed. This design allows data to be processed in a
single step, which is particularly useful for larger datasets. The
processing queue enables functions to be applied in a chunk-wise manner,
facilitating parallel processing and reducing memory demand.

Since the peaks data can be quite large, a processing queue is used to
ensure efficiency. Generally, the processing queue is applied either
temporarily when calling
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
or permanently when calling
[`applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html).
As explained below the processing efficiency can be further improved by
enabling chunk-wise processing.

## Usage

``` r
# S4 method for class 'Chromatograms'
applyProcessing(
  object,
  f = processingChunkFactor(object),
  BPPARAM = bpparam(),
  ...
)

# S4 method for class 'Chromatograms'
addProcessing(object, FUN, ...)

# S4 method for class 'Chromatograms'
processingChunkSize(object, ...)

# S4 method for class 'Chromatograms'
processingChunkSize(object) <- value

# S4 method for class 'Chromatograms'
processingChunkFactor(object, chunkSize = processingChunkSize(object), ...)
```

## Arguments

- object:

  A `Chromatograms` object.

- f:

  `factor` defining the grouping to split the `Chromatograms` object.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- ...:

  Additional arguments passed to the methods.

- FUN:

  For
  [`addProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html),
  a function to be added to the `Chromatograms` object's processing
  queue.

- value:

  `integer(1)` defining the chunk size.

- chunkSize:

  `integer(1)` for `processingChunkFactor` defining the chunk size. The
  default is the value stored in the `Chromatograms` object's
  `processingChunkSize` slot.

## Value

[`processingChunkSize()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
returns the currently defined processing chunk size (or `Inf` if it is
not defined).
[`processingChunkFactor()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
returns a `factor` defining the chunks into which `object` will be split
for (parallel) chunk-wise processing or a `factor` of length 0 if no
splitting is defined.

Refer to the individual function description for information on the
return value.

## Note

Some backends might not support parallel processing. For these, the
[`backendBpparam()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function will always return a
[`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
regardless of how parallel processing was defined.

## Apply Processing

The
[`applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
function applies the processing queue to the backend and returns the
updated `Chromatograms` object. The processing queue is a list of
processing steps applied to the chromatograms data. Each element in the
list is a function that processes the chromatograms data. To apply
processing to the peaks data, the backend must be set to a non-read-only
backend using the
[`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function.

## Parallel and Chunk-wise Processing of `Chromatograms`

Many operations on `Chromatograms` objects, especially those involving
the actual peaks data (see
[peaksData](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)),
support chunk-wise processing. This involves splitting the
`Chromatograms` into smaller parts (chunks) that are processed
iteratively. This enables parallel processing by data chunk and reduces
memory demand since only the peak data of the currently processed subset
is loaded into memory. Chunk-wise processing, which is disabled by
default, can be enabled by setting the processing chunk size of a
`Chromatograms` object using the
[`processingChunkSize()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
function to a value smaller than the length of the `Chromatograms`
object. For example, setting `processingChunkSize(chr) <- 1000` will
cause any data manipulation operation on `chr`, such as
[`filterPeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md),
to be performed in parallel for sets of 1000 chromatograms in each
iteration.

Chunk-wise processing is particularly useful for `Chromatograms` objects
using an *on-disk* backend or for very large experiments. For small
datasets or `Chromatograms` using an in-memory backend, direct
processing might be more efficient. Setting the chunk size to `Inf` will
disable chunk-wise processing.

Some backends may prefer a specific type of splitting and chunk-wise
processing. For example, the `ChromBackendMzR` backend needs to load MS
data from the original (mzML) files, so chunk-wise processing on a
per-file basis is ideal. The
[backendParallelFactor()](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
function for `ChromBackend` allows backends to suggest a preferred data
chunking by returning a `factor` defining the respective data chunks.
The `ChromBackendMzR` returns a `factor` based on the *dataOrigin*
chromatograms variable. A `factor` of length 0 is returned if no
particular preferred splitting is needed. The suggested chunk definition
will be used if no finite
[`processingChunkSize()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
is defined. Setting the `processingChunkSize` overrides
`backendParallelFactor`.

Functions to configure parallel or chunk-wise processing:

- [`processingChunkSize()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html):
  Gets or sets the size of the chunks for parallel or chunk-wise
  processing of a `Chromatograms` object. With a value of `Inf` (the
  default), no chunk-wise processing will be performed.

- [`processingChunkFactor()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html):
  Returns a `factor` defining the chunks into which a `Chromatograms`
  object will be split for chunk-wise (parallel) processing. A `factor`
  of length 0 indicates that no chunk-wise processing will be performed.

## Author

Johannes Rainer, Philippine Louail

## Examples

``` r
# Create a Chromatograms object
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    chromIndex = c(1L, 2L, 3L)
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

divide_intensities <- function(x, y, ...) {
    intensity(x) <- lapply(intensity(x), `/`, y)
    x
}

## Add the function to the procesing queue
chr <- addProcessing(chr, divide_intensities, y = 2)
chr
#> Chromatographic data (Chromatograms) with 3 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1          1       1 112.2
#> 2          2       1 123.3
#> 3          3       1 134.4
#> ... 0 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> Lazy evaluation queue: 1 processing step(s)

# Apply the processing queue
chr <- applyProcessing(chr)
```
