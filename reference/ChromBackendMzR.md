# Chromatographic Data Backend for Reading mzML Files

The `ChromBackendMzR` inherits all slots and methods from the base
`ChromBackendMemory` backend, providing additional functionality for
reading chromatographic data from mzML files.

Unlike the `ChromBackendMemory` backend, the `ChromBackendMzR` backend
should have the *dataOrigin* chromatographic variables populated with
the file path of the mzML file from which the chromatographic data was
read.

Note that the `ChromBackendMzR` backend is read-only and does not
support direct modification of chromatographic data. However, it does
support `peaksData` slot replacement, which will modify the `@peaksData`
slot but not the local mzML files. This is indicated by the "inMemory"
slot being set to TRUE.

Implementing functionalities with the `ChromBackendMzR` backend should
be simplified as much as possible and reuse the methods already
implemented for `ChromBackendMemory` when possible.

## Usage

``` r
ChromBackendMzR()

# S4 method for class 'ChromBackendMzR'
backendInitialize(object, files = character(), BPPARAM = bpparam(), ...)
```

## Arguments

- object:

  A `ChromBackendMzR` object.

- files:

  A character vector of file paths to mzML files.

- BPPARAM:

  Parallel setup configuration. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- ...:

  Additional parameters to be passed.

## Value

Refer to the individual function description for information on the
return value.

## Author

Philippine Louail

## Examples

``` r
library(mzR)
#> Loading required package: Rcpp
library(MsDataHub)

## Load an mzML file
MRM_file <- MRM.standmix.5.mzML()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

## Initialize the ChromBackendMzR object
be_empty <- ChromBackendMzR()
be <- backendInitialize(be_empty, files = MRM_file, BPPARAM = SerialParam())
```
