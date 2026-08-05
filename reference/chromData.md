# Chromatographic Peaks Metadata.

As explained in the
[`Chromatograms`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
class documentation, the `Chromatograms` object is a container for
chromatogram data that includes chromatographic peaks data (*retention
time* and related intensity values, also referred to as *peaks data
variables* in the context of `Chromatograms`) and metadata of individual
chromatograms (so called *chromatograms variables*).

The *chromatograms variables* information can be accessed using the
`chromData()` function. it is also possible to access specific
chromatograms variables using `$`.

`@chromData` can be accessed, replaced but also filtered/subsetted.
Refer to the sections below for more details.

## Usage

``` r
# S4 method for class 'Chromatograms'
chromData(object, columns = chromVariables(object), drop = FALSE)

# S4 method for class 'Chromatograms'
chromData(object) <- value

# S4 method for class 'Chromatograms'
chromVariables(object)

# S4 method for class 'Chromatograms'
chromIndex(object)

# S4 method for class 'Chromatograms'
chromIndex(object) <- value

# S4 method for class 'Chromatograms'
collisionEnergy(object)

# S4 method for class 'Chromatograms'
collisionEnergy(object) <- value

# S4 method for class 'Chromatograms'
dataOrigin(object)

# S4 method for class 'Chromatograms'
dataOrigin(object) <- value

# S4 method for class 'Chromatograms'
msLevel(object)

# S4 method for class 'Chromatograms'
msLevel(object) <- value

# S4 method for class 'Chromatograms'
mz(object)

# S4 method for class 'Chromatograms'
mz(object) <- value

# S4 method for class 'Chromatograms'
mzMax(object)

# S4 method for class 'Chromatograms'
mzMax(object) <- value

# S4 method for class 'Chromatograms'
mzMin(object)

# S4 method for class 'Chromatograms'
mzMin(object) <- value

# S4 method for class 'Chromatograms'
length(x)

# S4 method for class 'Chromatograms'
precursorMz(object)

# S4 method for class 'Chromatograms'
precursorMz(object) <- value

# S4 method for class 'Chromatograms'
precursorMzMin(object)

# S4 method for class 'Chromatograms'
precursorMzMin(object) <- value

# S4 method for class 'Chromatograms'
precursorMzMax(object)

# S4 method for class 'Chromatograms'
precursorMzMax(object) <- value

# S4 method for class 'Chromatograms'
productMz(object)

# S4 method for class 'Chromatograms'
productMz(object) <- value

# S4 method for class 'Chromatograms'
productMzMin(object)

# S4 method for class 'Chromatograms'
productMzMin(object) <- value

# S4 method for class 'Chromatograms'
productMzMax(object)

# S4 method for class 'Chromatograms'
productMzMax(object) <- value

# S4 method for class 'Chromatograms'
filterChromData(
  object,
  variables = character(),
  ranges = numeric(),
  match = c("any", "all"),
  keep = TRUE
)
```

## Arguments

- object:

  A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object.

- columns:

  A `character` vector of chromatograms variables to extract.

- drop:

  A `logical` indicating whether to drop dimensions when extracting a
  single variable.

- value:

  replacement value for `<-` methods. See individual method description
  or expected data type.

- x:

  A
  [Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  object.

- variables:

  For `filterChromData()`: `character` vector with the names of the
  chromatogram variables to filter for. The list of available
  chromatogram variables can be obtained with `chromVariables()`.

- ranges:

  For `filterChromData()` : a `numeric` vector of paired values (upper
  and lower boundary) that define the ranges to filter the `object`.
  These paired values need to be in the same order as the `variables`
  parameter (see below).

- match:

  For `filterChromData()` : `character(1) ` defining whether the
  condition has to match for all provided `ranges` (`match = "all"`; the
  default), or for any of them (`match = "any"`) for chromatogram data
  to be retained.

- keep:

  For `filterChromData()`: `logical(1)` defining whether to keep
  (`keep = TRUE`) or remove (`keep = FALSE`) the chromatogram data that
  match the condition.

## Value

Refer to the individual function description for information on the
return value.

## Chromatograms variables and accessor functions

The following chromatograms variables are guaranteed to be provided by a
`Chromatograms` object and to be accessible with either the
`chromData()` or a specific function named after the variables names:

- `chromIndex`: an `integer` with the index of the chromatogram in the
  original source file (e.g. *mzML* file).

- `collisionEnergy`: for SRM data, `numeric` with the collision energy
  of the precursor.

- `dataOrigin`: optional `character` with the origin of the data.

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

## Filter Chromatograms variables

Functions that filter `Chromatograms` based on chromatograms variables
(i.e, `@chromData` ) will remove chromatographic data that do not meet
the specified conditions. This means that if a chromatogram is filtered
out, its corresponding `@chromData` and `@peaksData` will be removed
from the object immediately.

The available functions to filter chromatogram data are:

- `filterChromData()`: Filters numerical chromatographic data variables
  based on the provided numerical `ranges`. The method returns a
  `Chromatograms` object containing only the chromatograms that match
  the specified conditions. This function results in an object with
  fewer chromatograms than the original.

## See also

[Chromatograms](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
for a general description of the `Chromatograms` object.
[peaksData](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
for a general description of the chromatographic peaks data available in
the object, as well as how to access, replace and subset them.
[processingQueue](https://rformassspectrometry.github.io/Chromatograms/reference/processingQueue.md)
for more information on the queuing of processings and parallelization
for larger dataset processing.

## Author

Philippine Louail

## Examples

``` r

# Create a Chromatograms object
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    chromIndex = c(1L, 2L, 3L)
)

be <- backendInitialize(new("ChromBackendMemory"), chromData = cdata)

chr <- Chromatograms(be)

# Access chromatograms variables
chromData(chr)
#>   msLevel    mz chromIndex collisionEnergy dataOrigin mzMin mzMax precursorMz
#> 1       1 112.2          1              NA       <NA>    NA    NA          NA
#> 2       1 123.3          2              NA       <NA>    NA    NA          NA
#> 3       1 134.4          3              NA       <NA>    NA    NA          NA
#>   precursorMzMin precursorMzMax productMz productMzMin productMzMax
#> 1             NA             NA        NA           NA           NA
#> 2             NA             NA        NA           NA           NA
#> 3             NA             NA        NA           NA           NA

# Access specific chromatograms variables
chromData(chr, columns = "msLevel")
#>   msLevel
#> 1       1
#> 2       1
#> 3       1

msLevel(chr)
#> [1] 1 1 1

# Replace chromatograms variables
msLevel(chr) <- c(1L, 2L, 2L)

# Filter chromatograms variables
filterChromData(chr,
    variables = "msLevel", ranges = c(1L, 1L),
    keep = FALSE
)
#> Chromatographic data (Chromatograms) with 2 chromatograms in a ChromBackendMemory backend:
#>   chromIndex msLevel    mz
#> 2          2       2 123.3
#> 3          3       2 134.4
#> ... 10 more  chromatogram variables/columns
#> ... 2 peaksData variables
```
