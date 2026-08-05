# Improved in-memory Chromatographic data backend

`ChromBackendMemory`: This backend stores chromatographic data directly
in memory, making it ideal for small datasets or testing. It can be
initialized with a `data.frame` of chromatographic data via the
`chromData` parameter and a `list` of `data.frame` entries for peaks
data using the `peaksData` parameter. These data can be accessed with
the
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
functions.

## Usage

``` r
ChromBackendMemory()

# S4 method for class 'ChromBackendMemory'
backendInitialize(
  object,
  chromData = fillCoreChromVariables(data.frame()),
  peaksData = list(.EMPTY_PEAKS_DATA),
  ...
)
```

## Arguments

- object:

  A `ChromBackendMemory` object.

- chromData:

  For `backendInitialize()` of a `ChromBackendMemory` backend, a
  `data.frame` with the chromatographic data. If not provided (or if
  empty), a default `data.frame` with the core chromatographicvariables
  will be created.

- peaksData:

  For `backendInitialize()` of a `ChromBackendMemory` backend, a `list`
  of `data.frame` with the peaks data. If not provided (or if empty), a
  default `list` of empty `data.frame` with the core peaks variables
  will be created. The length of the list should match the number of
  chromatograms in the `chromData` parameter.

- ...:

  Additional parameters to be passed.

## Value

Refer to the individual function description for information on the
return value.

## Author

Philippine Louail

## Examples

``` r

## Method 1: Initialize backend directly
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

cbm <- ChromBackendMemory()
cbm <- backendInitialize(cbm, chromData = cdata, peaksData = pdata)
cbm
#> ChromBackendMemory with 3 chromatograms
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       1 134.4
#> ... 3 more  chromatogram variables/columns
#> ... 2 peaksData variables

## Method 2: Use Chromatograms constructor (recommended)
chr <- Chromatograms(ChromBackendMemory(), chromData = cdata, peaksData = pdata)
chr
#> Chromatographic data (Chromatograms) with 3 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 1         NA       1 112.2
#> 2         NA       1 123.3
#> 3         NA       1 134.4
#> ... 3 more  chromatogram variables/columns
#> ... 2 peaksData variables
```
