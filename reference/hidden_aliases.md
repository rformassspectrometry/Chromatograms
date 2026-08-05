# Fill data.frame with columns for missing core chromatogram variables.

`fillCoreChromVariables()` fills a provided `data.frame` with columns
for eventually missing *core* chromatogram variables. The missing core
variables are added as new columns with missing values (`NA`) of the
correct data type. Use
[`coreChromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
to list the set of core variables and their data types.

`validChromData()` checks that columns, representing *core* chromatogram
variables are of the correct data type.

For S4 methods that require a documentation entry but only clutter the
index.

This method returns the chromatographic data stored in the backend. If
not specified otherwise it will return all defined columns in the
chromData slot as well as adding the `coreChromVariables` missing with
NA values.

## Usage

``` r
reset(object, ...)

plotChromatogramsOverlay(object, ...)

fillCoreChromVariables(x = data.frame())

validChromData(x = data.frame(), error = TRUE)

validPeaksData(x = list(), error = TRUE)

# S4 method for class 'ChromBackendMemory'
backendMerge(object, ...)

# S4 method for class 'ChromBackendMemory'
chromData(object, columns = chromVariables(object), drop = FALSE)

# S4 method for class 'ChromBackendMemory'
chromData(object) <- value

# S4 method for class 'ChromBackendMemory'
chromVariables(object)

# S4 method for class 'ChromBackendMemory'
peaksData(object, columns = peaksVariables(object), drop = FALSE, ...)

# S4 method for class 'ChromBackendMemory'
peaksData(object) <- value

# S4 method for class 'ChromBackendMemory'
peaksVariables(object)

# S4 method for class 'ChromBackendMemory'
backendParallelFactor(object, ...)

# S4 method for class 'ChromBackendMemory'
isReadOnly(object)

# S4 method for class 'ChromBackendMemory'
lengths(x)

# S4 method for class 'ChromBackendMemory'
intensity(object)

# S4 method for class 'ChromBackendMemory'
rtime(object)

# S4 method for class 'ChromBackendMemory'
show(object)

# S4 method for class 'ChromBackendMemory'
supportsSetBackend(object, ...)

# S4 method for class 'ChromBackendMemory'
x[i, j, ..., drop = FALSE]

# S4 method for class 'ChromBackendMemory'
x$name

# S4 method for class 'ChromBackendMemory'
x$name <- value

# S4 method for class 'ChromBackendMemory'
chromExtract(object, peak.table, by)

# S4 method for class 'ChromBackendMzR'
show(object)

# S4 method for class 'ChromBackendMzR'
backendParallelFactor(object, ...)

# S4 method for class 'ChromBackendMzR'
isReadOnly(object)

# S4 method for class 'ChromBackendMzR'
peaksData(
  object,
  columns = peaksVariables(object),
  drop = FALSE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'ChromBackendMzR'
peaksData(object) <- value

# S4 method for class 'ChromBackendMzR'
chromData(object) <- value

# S4 method for class 'ChromBackendMzR'
supportsSetBackend(object, ...)

# S4 method for class 'ChromBackendMzR'
intensity(object)

# S4 method for class 'ChromBackendMzR'
rtime(object)

# S4 method for class 'ChromBackendMzR'
lengths(x)

# S4 method for class 'ChromBackendMzR'
x[i, j, ..., drop = TRUE]

# S4 method for class 'ChromBackendMzR'
chromExtract(object, peak.table, by, ...)

# S4 method for class 'ChromBackendSpectra'
show(object)

# S4 method for class 'ChromBackendSpectra'
factorize(object, factorize.by = c("msLevel", "dataOrigin"), ...)

# S4 method for class 'ChromBackendSpectra'
backendParallelFactor(object, ...)

# S4 method for class 'ChromBackendSpectra'
isReadOnly(object)

# S4 method for class 'ChromBackendSpectra'
peaksData(object, columns = peaksVariables(object), drop = FALSE, ...)

# S4 method for class 'ChromBackendSpectra'
peaksData(object) <- value

# S4 method for class 'ChromBackendSpectra'
supportsSetBackend(object, ...)

# S4 method for class 'ChromBackendSpectra'
intensity(object)

# S4 method for class 'ChromBackendSpectra'
rtime(object)

# S4 method for class 'ChromBackendSpectra'
lengths(x)

# S4 method for class 'ChromBackendSpectra'
x[i, j, ..., drop = TRUE]

# S4 method for class 'ChromBackendSpectra'
chromExtract(object, peak.table, by, ...)

# S4 method for class 'Chromatograms'
show(object)
```

## Arguments

- object:

  A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object.

- x:

  `list` representing the peaks data of a `Chromatograms`

- error:

  `logical(1)` whether an error should be thrown (the default) if one or
  more columns don't have the correct data type.

## Value

input data frame `x` with missing core variables added (with the correct
data type).

If the core variables have all the correct data type: an empty
character. If one or more core variables (columns) have the wrong data
type the function either throws an error (with `error = TRUE`) or
returns a `character` specifying which variables/columns don't have the
correct type (for `error = FALSE`).

Not applicable

## Examples

``` r

## Define a data frame
a <- data.frame(msLevel = c(1L, 1L, 2L), other_column = "b")

## Add missing core chromatogram variables to this data frame
fillCoreChromVariables(a)
#>   msLevel other_column chromIndex collisionEnergy dataOrigin mz mzMin mzMax
#> 1       1            b         NA              NA       <NA> NA    NA    NA
#> 2       1            b         NA              NA       <NA> NA    NA    NA
#> 3       2            b         NA              NA       <NA> NA    NA    NA
#>   precursorMz precursorMzMin precursorMzMax productMz productMzMin productMzMax
#> 1          NA             NA             NA        NA           NA           NA
#> 2          NA             NA             NA        NA           NA           NA
#> 3          NA             NA             NA        NA           NA           NA

## The data.frame thus contains columns for all core chromatogram
## variables in the respective expected data type (but filled with
## missing values).
```
