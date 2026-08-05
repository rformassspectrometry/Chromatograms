# Merging, combining and splitting Chromatograms

Various functions are available to combine or split data from one or
more `Chromatograms` objects. These are:

- [`c()`](https://rdrr.io/r/base/c.html) and
  `concatenateChromatograms()`: combines several `Chromatograms` objects
  into a single object. The resulting `Chromatograms` contains all data
  from all individual `Chromatograms`, i.e. the union of all their
  chromatograms variables. Concatenation will fail if the processing
  queue of any of the `Chromatograms` objects is not empty or if
  different backends are used for the `Chromatograms` objects. In such
  cases it is suggested to first change the backends of all
  `Chromatograms` to the same type of backend (using the
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  function) and to eventually (if needed) apply the processing queue
  using the
  [`applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
  function.

- [`split()`](https://rdrr.io/r/base/split.html): splits the
  `Chromatograms` object based on a provided grouping factor returning a
  `list` of `Chromatograms` objects.

## Usage

``` r
concatenateChromatograms(x, ...)

# S4 method for class 'Chromatograms'
c(x, ...)

# S4 method for class 'Chromatograms,ANY'
split(x, f, drop = FALSE, ...)
```

## Arguments

- x:

  A `Chromatograms` object.

- ...:

  For [`c()`](https://rdrr.io/r/base/c.html) and
  `concatenateChromatograms()`: `Chromatograms` objects or a `list` of
  `Chromatograms` objects.

- f:

  `factor` defining the grouping to split the `Chromatograms` object.

- drop:

  For [`split()`](https://rdrr.io/r/base/split.html): not considered.

## Value

- [`c()`](https://rdrr.io/r/base/c.html) and
  `concatenateChromatograms()`: a single `Chromatograms` object
  containing the data from all input objects.

- [`split()`](https://rdrr.io/r/base/split.html): a `list` of
  `Chromatograms` objects.

## Examples

``` r

## Create two Chromatograms objects
cdata1 <- data.frame(
    msLevel = c(1L, 1L),
    mz = c(112.2, 123.3),
    dataOrigin = c("file1", "file1")
)
pdata1 <- list(
    data.frame(rtime = c(1.0, 2.0, 3.0), intensity = c(100, 200, 150)),
    data.frame(rtime = c(1.0, 2.0, 3.0), intensity = c(80, 120, 90))
)
chr1 <- Chromatograms(
    ChromBackendMemory(),
    chromData = cdata1,
    peaksData = pdata1
)

cdata2 <- data.frame(
    msLevel = c(2L, 2L),
    mz = c(134.4, 145.5),
    dataOrigin = c("file2", "file2")
)
pdata2 <- list(
    data.frame(rtime = c(4.0, 5.0, 6.0), intensity = c(300, 400, 350)),
    data.frame(rtime = c(4.0, 5.0, 6.0), intensity = c(200, 250, 180))
)
chr2 <- Chromatograms(
    ChromBackendMemory(),
    chromData = cdata2,
    peaksData = pdata2
)

## Combine using c()
chr_combined <- c(chr1, chr2)
chr_combined
#> Chromatographic data (Chromatograms) with 4 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       2 134.4
#> 4         NA       2 145.5
#> ... 4 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> Processing:
#>  Merged 2 Chromatograms into one [Wed Aug  5 11:48:05 2026] 

## Combine using concatenateChromatograms
chr_combined2 <- concatenateChromatograms(chr1, chr2)

## Combine a list of Chromatograms
chr_list <- list(chr1, chr2)
chr_combined3 <- concatenateChromatograms(chr_list)

## Split by msLevel
chr_split <- split(chr_combined, f = chr_combined$msLevel)
chr_split
#> $`1`
#> Chromatographic data (Chromatograms) with 2 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> ... 2 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> Processing:
#>  Merged 2 Chromatograms into one [Wed Aug  5 11:48:05 2026] 
#> 
#> $`2`
#> Chromatographic data (Chromatograms) with 2 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 3         NA       2 134.4
#> 4         NA       2 145.5
#> ... 2 more  chromatogram variables/columns
#> ... 2 peaksData variables
#> Processing:
#>  Merged 2 Chromatograms into one [Wed Aug  5 11:48:05 2026] 
#> 
```
