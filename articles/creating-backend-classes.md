# Creating new \`ChromBackend\` classes for Chromatograms

**Package**: Chromatograms 1.3.3\
**Compiled**: Wed Aug 5 11:48:09 2026

## Introduction

Similar to the
*[Spectra](https://bioconductor.org/packages/3.24/Spectra)* package, the
*[Chromatograms](https://bioconductor.org/packages/3.24/Chromatograms)*
also separates the user-faced functionality to process and analyze
chromatographic mass spectrometry (MS) data from the code for storage
and *representation* of the data. The latter functionality is provided
by implementations of the `ChromBackend` class, further on called
*backends*. This vignette describes the `ChromBackend` class and
illustrates on a simple example how a backend extending this class could
be implemented.

Contributions to this vignette (content or correction of typos) or
requests for additional details and information are highly welcome,
ideally *via* pull requests or *issues* on the package’s [github
repository](https://github.com/RforMassSpectrometry/Chromatograms).

This vignette describe the structure of a backend class and the methods
that need to be implemented. In order to see the structure of the
backend, the `@` accessor is used to access the slots of the backend
object. This is possible because the backend class is defined as an S4
class. However users should not use the `@` accessor to access the slots
of a backend object, but instead use the methods defined by the
`ChromBackend` class.

## What is a `ChromBackend`?

The purpose of a backend class extending the virtual `ChromBackend` is
to provide the chromatographic MS data to the `Chromatograms` object,
which is used by the user to interact with - and analyze the data. The
`ChromBackend` defines the API that new backends need to provide so that
they can be used with `Chromatograms`. This API defines a set of methods
to access the data. For many functions default implementations exist and
a dedicated implementation for a new backend is only needed if necessary
(e.g. if the data is stored in a way that a different access to it would
be better). In addition, a core set of variables (data fields), the so
called *core* chromatogram variables, is defined to describe the
chromatographic data. Each backend needs to provide these, but can also
define additional data fields. Before implementing a new backend it is
highly suggested to carefully read the following *Conventions and
definitions* section.

### Conventions and definitions

General conventions for chromatographic MS data of a `Chromatograms`
are:

- One `Chromatograms` object is designed to contain multiple
  chromatographic data (not data from a single chromatogram).
- retention time values within each chromatogram are expected to be
  sorted increasingly.
- Missing values (`NA`) for retention time values are not supported.
- Properties (data fields) of a chromatogram are called *chromatogram
  variables*. While backends can define their own properties, a minimum
  required set of chromatogram variables **must** be provided by each
  backend (even if their values are empty). These *core chromatogram
  variables* are listed (along with their expected data type) by the
  [`coreChromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
  function.
- `dataOrigin` defines for each chromatogram where the data is from
  expected to be of type`character`. Missing values should be
  `NA_character_`
- `ChromBackend` implementations can also represent purely *read-only*
  data resources. In this case only data accessor methods need to be
  implemented but not data replacement methods (i.e. `<-` methods that
  would allow to add or set variables. Read-only backends should
  implement the
  [`isReadOnly()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  method, that should then return `TRUE`.

### Notes on parallel and chunk-wise processing

For parallel processing, `Chromatograms` splits the backend based on a
defined `factor` and processes each in parallel (or *in serial* if a
`SerialParam` is used). The splitting `factor` can be defined for
`Chromatograms` by setting the parameter `processingChunkSize`.
Alternatively, through the
[`backendParallelFactor()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method the backend can also *suggest* a `factor` that should/could be
used for splitting and parallel processing. The default implementation
for
[`backendParallelFactor()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
is to return an empty `factor`
([`factor()`](https://rdrr.io/r/base/factor.html)) hence not suggesting
any preferred splitting. Besides parallel processing, for on-disk
backends (i.e., backends that don’t keep all of the data in memory),
this chunk-wise processing can also reduce the memory demand for
operations, because only the peak data of the current chunk needs to be
realized in memory.

## API

The `ChromBackend` class defines core methods that have to be
implemented by a MS *backend* as well as *optional* methods for which a
default implementation is already available. These functions are
described in sections *Required methods* and *Optional methods*,
respectively.

To create a new backend a class extending the virtual `ChromBackend`
needs to be implemented. In the following example we define a simple
class that uses a `data.frame` to store general properties
(*chromatogram variables*) and a list of `data.frame` for the retention
time and intensity values of each chromatograms, which represent the
actual chromatographic MS data. These values are store in a`list`, where
each element correspond to one chromatogram, as the number of values
(*peaks*) can vary between chromatograms. We also provide a basic
constructor function that returns an empty instance of the new class.

``` r

library(Chromatograms)

#' Definition of the backend class extending ChromBackend
setClass("ChromBackendTest",
    contains = "ChromBackend",
    slots = c(
        chromData = "data.frame",
        peaksData = "list"
    ),
    prototype = prototype(
        chromData = data.frame(),
        peaksData = list()
    )
)

#' Simple constructor function
ChromBackendTest <- function() {
    new("ChromBackendTest")
}
```

The 2 slots `@chromData` and `@peaksData` will be used to store the
general properties of the chromatograms and the actual chromatographic
data, respectively. each row in `chromData` will contain data for one
chromatogram with the columns being the different *chromatogram
variables* (i.e. additional properties of a chromatogram such as its m/z
value or MS level) and each element in `@peaksData` a `data.frame` with
the retention time and intensity values representing thus the *peaks*
data of the respective chromatogram. This is only one of the possibly
many ways chromatographic data might be represented.

We should ideally also add some basic validity function that ensures the
data to be correct (valid). The function below simply checks that the
number of rows of the `@chromData` slot matches the length of the
`@peaksData` slots.

``` r

#' Basic validation function
setValidity("ChromBackendTest", function(object) {
    if (length(object@peaksData) != nrow(object@chromData)) {
        return(
            "length of 'peaksData' has to match the number of rows of ",
            "'chromData'"
        )
    }
    NULL
})
```

    ## Class "ChromBackendTest" [in ".GlobalEnv"]
    ## 
    ## Slots:
    ##                                        
    ## Name:   chromData  peaksData    version
    ## Class: data.frame       list  character
    ## 
    ## Extends: 
    ## Class "ChromBackend", directly
    ## Class "ChromBackendOrMissing", by class "ChromBackend", distance 2

We can now create an instance of our new class with the
`ChromBackendTest()` function.

``` r

#' Create an empty instance of ChromBackendTest
be <- ChromBackendTest()
be
```

    ## An object of class "ChromBackendTest"
    ## Slot "chromData":
    ## data frame with 0 columns and 0 rows
    ## 
    ## Slot "peaksData":
    ## list()
    ## 
    ## Slot "version":
    ## [1] "0.1"

A [`show()`](https://rdrr.io/r/methods/show.html) method would allow for
a more convenient way how general information of our object is
displayed. Below we add an implementation of the
[`show()`](https://rdrr.io/r/methods/show.html) method.

``` r

#' implementation of show for ChromBackendTest
setMethod("show", "ChromBackendTest", function(object) {
    cd <- object@chromData
    cat(class(object), "with", nrow(cd), "chromatograms\n")
})
be
```

    ## ChromBackendTest with 0 chromatograms

### Required methods

Methods listed in this section **must** be implemented for a new class
extending `ChromBackend`. Methods should ideally also be implemented in
the order they are listed here. Also, it is strongly advised to write
dedicated unit tests for each newly implemented method or function
already **during** the development.

#### `dataStorage()`

The `dataStorage` chromatogram variable provides information how or
where the data is stored. The
[`dataStorage()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
method should therefore return a `character` vector of length equal to
the number of chromatograms that are represented by the object. The
values for `dataStorage` can be any character value, except `NA`. For
our example backend we define a simple
[`dataStorage()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
method that simply returns the column `"dataStorage"` from the
`@chromData` (as a `character`).

``` r

#' dataStorage method to provide information *where* data is stored
setMethod("dataStorage", "ChromBackendTest", function(object) {
    as.character(object@chromData$dataStorage)
})
```

Calling
[`dataStorage()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
on our example backend will thus return an empty `character` (since the
object created above does not contain any data).

``` r

dataStorage(be)
```

    ## character(0)

#### `length()`

[`length()`](https://rdrr.io/r/base/length.html) is expected to return
an `integer` of length 1 with the total number of chromatograms that are
represented by the backend. For our example backend we simply return the
number of rows of the `data.frame` stored in the `@chromData` slot.

``` r

#' length to provide information on the number of chromatograms
setMethod("length", "ChromBackendTest", function(x) {
    nrow(x@chromData)
})
length(be)
```

    ## [1] 0

#### `backendInitialize()`

The
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method should be called after creating an instance of the backend class
and is responsible for preparing (initializing) the backend with data.
This method can accept any parameters required by the backend to load or
initialize the data, such as file names, a database connection, or
objects containing the data. It is also recommended that the the special
chromatogram variables `dataStorage` and `dataOrigin` are set during
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html).

It is strongly recommended to validate the input data within the
initialize method. The advantage of performing these validity checks in
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
rather than using
[`setValidity()`](https://rdrr.io/r/methods/validObject.html) is that
computationally expensive operations/checks would only be performed
once,during initialization, instead of each time values within the
object are modified (e.g., through subsetting or similar operations),
which would occur with
[`setValidity()`](https://rdrr.io/r/methods/validObject.html).

We also use the
[`validChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
and
[`validPeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
functions to ensure that core chromatogram variables and core peaks
variables have the correct data type. These checks verify that
the`peaksData` contains only numeric values and that the number of
retention time and intensity values matches for each chromatogram.

Below we define a
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method that accepts a `data.frame` containing chromatogram variables and
a `list` with retention time and intensity values for each chromatogram.

``` r

#' backendInitialize method to fill the backend with data.
setMethod(
    "backendInitialize", "ChromBackendTest",
    function(object, chromData, peaksData) {
        if (!is.data.frame(chromData)) {
            stop(
                "'chromData' needs to be a 'data.frame' with the general",
                "chromatogram variables"
            )
        }
        ## Defining dataStorage and dataOrigin, if not available
        if (is.null(chromData$dataOrigin)) {
            chromData$dataOrigin <- NA_character_
        }
        ## Validate the provided data
        validChromData(chromData)
        validPeaksData(peaksData)
        ## Fill the object with data
        object@chromData <- chromData
        object@peaksData <- peaksData
        object
    }
)
```

In addition to adding the data to object, the function also define the
`dataOrigin` chromatographic variables. This variable is expected to
provide information on where the data is originating.

We can now create an instance of our backend class and fill it with
data. We thus first define our MS data and pass this to the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method.

``` r

# A data.frame with chromatogram variables.
cdata <- data.frame(
    msLevel = c(1L, 1L),
    mz = c(112.2, 123.3)
)

# Retention time and intensity values for each chromatogram.
pdata <- list(
    data.frame(
        rtime = c(12.4, 12.8, 13.2, 14.6),
        intensity = c(123.3, 153.6, 2354.3, 243.4)
    ),
    data.frame(
        rtime = c(45.1, 46.2),
        intensity = c(100, 80.1)
    )
)

#' Create and initialize the backend
be <- backendInitialize(ChromBackendTest(),
    chromData = cdata, peaksData = pdata
)
be
```

    ## ChromBackendTest with 2 chromatograms

This
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
implementation should assure data validity and integrity. Below we use
this function again to create our backend instance.

The
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method that we implemented for our backend class expects the user to
provide the full MS data. It would alternatively also be possible to
implement a method that takes data file names as input from which the
function can then import the data. The purpose of the
[`backendInitialize()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method is to *initialize* and prepare the data in a way that it can be
accessed by a `Chromatograms` object. Whether the data is actually
loaded into memory or simply referenced and loaded upon request does not
matter as long as the backend is able to provide the data though its
accessor methods when requested by the `Chromatograms` object.

#### `chromVariables()`

The
[`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
method should return a `character` vector with the names of all
available chromatogram variables of the backend. While a backend class
should support defining and providing their own variables, each
`ChromBackend` class **must** provide also the *core chromatogram
variables* (in the correct data type). These can be listed by the
[`coreChromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
function:

``` r

#' List core chromatogram variables along with data types.
coreChromVariables()
```

    ##      chromIndex collisionEnergy      dataOrigin         msLevel              mz 
    ##       "integer"       "numeric"     "character"       "integer"       "numeric" 
    ##           mzMin           mzMax     precursorMz  precursorMzMin  precursorMzMax 
    ##       "numeric"       "numeric"       "numeric"       "numeric"       "numeric" 
    ##       productMz    productMzMin    productMzMax 
    ##       "numeric"       "numeric"       "numeric"

A typical
[`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
method for a `ChromBackend` class will thus be implemented similarly to
the one for our `ChromBackendTest` test backend: it will return the
names for all available chromatogram variables that can be called by
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
within the backend object. There is a default implementation for
[`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
that will return the core chromatogram variables. However if a backend
class defines additional chromatogram variables, the
[`chromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
method should be implemented to return the names of these additional
variables as well.

``` r

#' Accessor for available chromatogram variables
setMethod("chromVariables", "ChromBackendTest", function(object) {
    union(names(object@chromData), names(coreChromVariables()))
})

chromVariables(be)
```

    ##  [1] "msLevel"         "mz"              "dataOrigin"      "chromIndex"     
    ##  [5] "collisionEnergy" "mzMin"           "mzMax"           "precursorMz"    
    ##  [9] "precursorMzMin"  "precursorMzMax"  "productMz"       "productMzMin"   
    ## [13] "productMzMax"

#### `chromData()`

The `chromData` method should return the **full** chromatogram data
within a backend as a `data.frame` object. A parameter `columns` should
allow to define the names of the variables that should be returned. A
parameter `drop` should also be implemented to allow for the calling of
one column while still controlling the return type. Each row in this
data frame should represent one chromatogram, each column a chromatogram
variable. The `data.frame` **must** provide values (even if they are
`NA`) for **all** requested chromatogram variables of the backend
(**including** the core chromatogram variables). The
[`fillCoreChromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
function from the *Chromatograms* package allows to *complete* (fill) a
provided `data.frame` with eventually missing core chromatogram
variables:

``` r

#' Get the data.frame with the available chrom variables
be@chromData
```

    ##   msLevel    mz dataOrigin
    ## 1       1 112.2       <NA>
    ## 2       1 123.3       <NA>

``` r

#' Complete this data.frame with missing core variables
fillCoreChromVariables(be@chromData)
```

    ##   msLevel    mz dataOrigin chromIndex collisionEnergy mzMin mzMax precursorMz
    ## 1       1 112.2       <NA>         NA              NA    NA    NA          NA
    ## 2       1 123.3       <NA>         NA              NA    NA    NA          NA
    ##   precursorMzMin precursorMzMax productMz productMzMin productMzMax
    ## 1             NA             NA        NA           NA           NA
    ## 2             NA             NA        NA           NA           NA

We can thus use this function to add eventually missing core
chromatogram variables in the `chromData` implementation for our
backend:

``` r

#' function to extract the full chromData
setMethod(
    "chromData", "ChromBackendTest",
    function(object, columns = chromVariables(object),
    drop = FALSE) {
        if (!any(chromVariables(object) %in% columns)) {
            stop(
                "Some of the requested Chromatogram variables are not ",
                "available"
            )
        }
        res <- fillCoreChromVariables(object@chromData)
        res <- res[, columns, drop = drop]
        res
    }
)
```

We can now use
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
to either extract the full chromatogram data from the backend, or only
the data for selected variables.

``` r

#' Extract the full data
chromData(be)
```

    ##   msLevel    mz dataOrigin chromIndex collisionEnergy mzMin mzMax precursorMz
    ## 1       1 112.2       <NA>         NA              NA    NA    NA          NA
    ## 2       1 123.3       <NA>         NA              NA    NA    NA          NA
    ##   precursorMzMin precursorMzMax productMz productMzMin productMzMax
    ## 1             NA             NA        NA           NA           NA
    ## 2             NA             NA        NA           NA           NA

``` r

#' Selected variables
chromData(be, c("mz", "msLevel"))
```

    ##      mz msLevel
    ## 1 112.2       1
    ## 2 123.3       1

``` r

#' Only missing core chromatograms variables
chromData(be, c("collisionEnergy", "mzMin"))
```

    ##   collisionEnergy mzMin
    ## 1              NA    NA
    ## 2              NA    NA

#### `peaksVariables()`

The
[`peaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
function is supposed to provide the names of the available *peaks
variables*. If additional peaks variables would be available, these
could also be listed by the
[`peaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method. There is a default implementation for `peaksVaraibles()` that
will return the core peaks variables. However if a backend class defines
additional peaks variables, the
[`peaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method should be implemented to return the names of these additional
variables as well.

``` r

setMethod("peaksVariables", "ChromBackendTest", function(object) {
    union(names(corePeaksVariables()), names(object@peaksData[[1]]))
})
```

We can now see what peaks variables are present in our object:

``` r

peaksVariables(be)
```

    ## [1] "rtime"     "intensity"

#### `peaksData()`

The
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method extracts the chromatographic data (*peaks*), i.e., the
chromatograms’ retention time and intensity values. This data is
returned as a `list` of `data.frame`, with one array per chromatogram
with columns being the *peaks variables* (retention time and intensity
values) and rows the individual data pairs. Each backend must provide
retention times and intensity values with this method, but additional
peaks variables (columns) are also supported.

In a similar way as for the chromatogram variables, a backend should
support defining and providing their own variables and each
`ChromBackend` class **must** provide also the *core peaks variables*
(in the correct data type). These can be listed by the
[`corePeaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
function:

``` r

corePeaksVariables()
```

    ##     rtime intensity 
    ## "numeric" "numeric"

Below we implement the
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method for our backend.

``` r

#' method to extract the full chromatographic data as list of arrays
setMethod(
    "peaksData", "ChromBackendTest",
    function(object, columns = peaksVariables(object), drop = FALSE) {
        if (!all(columns %in% peaksVariables(object))) {
            stop("Some of the requested peaks variables are not available")
        }
        res <- lapply(object@peaksData, function(x) x[, columns, drop = drop])
        res
    }
)
```

And with this method we can now extract the peaks data from our backend.

``` r

#' Extract the *peaks* data (i.e. intensity and retention times)
peaksData(be)
```

    ## [[1]]
    ##   rtime intensity
    ## 1  12.4     123.3
    ## 2  12.8     153.6
    ## 3  13.2    2354.3
    ## 4  14.6     243.4
    ## 
    ## [[2]]
    ##   rtime intensity
    ## 1  45.1     100.0
    ## 2  46.2      80.1

Since the
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method is the main function used by a `Chromatograms` to retrieve data
from the backend (and further process the values), this method should be
implemented in an efficient way.

#### `[`

The `[` method allows to subset `ChromBackend` objects. This operation
is expected to reduce a `ChromBackend` object to the selected
chromatograms without changing values for the subset chromatograms. The
method should support to subset by indices or logical vectors and should
also support duplicating elements (i.e., when duplicated indices are
used) as well as to subset in arbitrary order. An error should be thrown
if indices are out of bounds, but the method should also support
returning an empty backend with `[integer()]`. The
[`MsCoreUtils::i2index`](https://rdrr.io/pkg/MsCoreUtils/man/i2index.html)
function can be used to check and convert the provided parameter `i`
(defining the subset) to an integer vector.

Below we implement a possible `[` for our test backend class. We ignore
the parameters `j` from the definition of the `[` generic, since we
treat our data to be one-dimensional (with each chromatogram being one
element).

``` r

#' Main subset method.
setMethod("[", "ChromBackendTest", function(x, i, j, ..., drop = FALSE) {
    i <- MsCoreUtils::i2index(i, length = length(x))
    x@chromData <- x@chromData[i, ]
    x@peaksData <- x@peaksData[i]
    x
})
```

We can now subset our backend to the last two chromatograms.

``` r

a <- be[1]
chromData(a)
```

    ##   msLevel    mz dataOrigin chromIndex collisionEnergy mzMin mzMax precursorMz
    ## 1       1 112.2       <NA>         NA              NA    NA    NA          NA
    ##   precursorMzMin precursorMzMax productMz productMzMin productMzMax
    ## 1             NA             NA        NA           NA           NA

Or extracting the second chromatogram multiple times.

``` r

a <- be[c(1, 1, 1)]
chromData(a)
```

    ##     msLevel    mz dataOrigin chromIndex collisionEnergy mzMin mzMax precursorMz
    ## 1         1 112.2       <NA>         NA              NA    NA    NA          NA
    ## 1.1       1 112.2       <NA>         NA              NA    NA    NA          NA
    ## 1.2       1 112.2       <NA>         NA              NA    NA    NA          NA
    ##     precursorMzMin precursorMzMax productMz productMzMin productMzMax
    ## 1               NA             NA        NA           NA           NA
    ## 1.1             NA             NA        NA           NA           NA
    ## 1.2             NA             NA        NA           NA           NA

#### `$`

The `$` method is expected to extract a single chromatogram or peaks
variable from a backend. Parameter `name` should allow to name the
variable to return. Each `ChromBackend` **must** support extracting the
core chromatogram and core peaks variables with this method (even if no
data might be available for that variable). In our example
implementation below we make use of the
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
method, but more efficient implementations might be possible as well.
Also, the `$` method should check if the requested variable is available
and should throw an error otherwise.

``` r

#' Access a single chromatogram variable
setMethod("$", "ChromBackendTest", function(x, name) {
    if (name %in% union(chromVariables(x), names(coreChromVariables()))) {
        res <- chromData(x, columns = name, drop = TRUE)
    } else if (name %in% peaksVariables(x)) {
        res <- peaksData(x, columns = name, drop = TRUE)
    } else {
        stop("The requested variable '", name, "' is not available")
    }
    res
})
```

With this we can now extract the MS levels

``` r

be$msLevel
```

    ## [1] 1 1

or a core chromatogram variable without values in our example backend.

``` r

be$precursorMz
```

    ## [1] NA NA

or also the intensity values

``` r

be$intensity
```

    ## [[1]]
    ## [1]  123.3  153.6 2354.3  243.4
    ## 
    ## [[2]]
    ## [1] 100.0  80.1

#### `backendMerge()`

The
[`backendMerge()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method merges (combines) `ChromBackend` objects (of the same type!) into
a single instance. For our test backend we thus need to combine the
values in the `@chromData`, `@peaksData` slots. To support also merging
of `data.frame`s with different sets of columns we use the
[`MsCoreUtils::rbindFill`](https://rdrr.io/pkg/MsCoreUtils/man/rbindFill.html)
function instead of a simple `rbind` (this function joins data frames
making an union of all available columns filling eventually missing
columns with `NA`).

``` r

#' Method allowing to join (concatenate) backends
setMethod("backendMerge", "ChromBackendTest", function(object, ...) {
    res <- object
    object <- unname(c(list(object), list(...)))
    res@peaksData <- do.call(c, lapply(object, function(z) z@peaksData))
    res@chromData <- do.call(
        MsCoreUtils::rbindFill,
        lapply(object, function(z) z@chromData)
    )
    validObject(res)
    res
})
```

Testing the function by merging the example backend instance with
itself.

``` r

a <- backendMerge(be, be[2], be)
a
```

    ## ChromBackendTest with 5 chromatograms

### Data replacement methods

As stated in the general description, `ChromBackend` implementations can
also be purely *read-only* resources allowing to just access, but not to
replace data. For these backends
[`isReadOnly()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
should return `FALSE`. Data replacement methods listed in this section
would not need to be implemented. Our example backend stores the full
data in memory, within the object, and hence we can easily change and
replace values.

Since we support replacing values we also implement the
[`isReadOnly()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method for our example implementation to return `FALSE` (instead of the
default `TRUE`).

``` r

#' Default for backends:
isReadOnly(be)
```

    ## [1] FALSE

``` r

#' Implementation of isReadOnly for ChromBackendTest
setMethod("isReadOnly", "ChromBackendTest", function(object) FALSE)
isReadOnly(be)
```

    ## [1] FALSE

All data replacement function are expected to return an instance of the
same backend class that was used as input.

#### `chromData<-`

The main replacement method is `chromData<-` which should allow to
replace the chormtaogram variables content of a backend with new data.
This data is expected to be provided as a `data.frame` (similar to the
one returned by
[`chromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)).
While values can be replaced, the number of chromatograms before and
after a call to `chromData<-` has to be the same.

``` r

#' Replacement method for the full chromatogram data
setReplaceMethod("chromData", "ChromBackendTest", function(object, value) {
    if (is(value, "DataFrame")) {
        value <- as(value, "data.frame")
    }
    if (!inherits(value, "data.frame")) {
        stop("'value' is expected to be a 'data.frame'")
    }
    if (length(object) && length(object) != nrow(value)) {
        stop("'value' has to be a 'data.frame' with ", length(object), " rows")
    }
    validChromData(value)
    object@chromData <- value
    object
})
```

To test this new method we extract the full chromatogram data from our
example data set, add an additional column (chromatogram variable) and
use `chromData<-` to replace the data of the backend.

``` r

d <- chromData(be)
d$new_col <- c("a", "b")

chromData(be) <- d
```

Check that we have now also the new column available.

``` r

be$new_col
```

    ## [1] "a" "b"

#### `$<-`

The `$<-` method should allow to replace values for an existing
chromatogram variable or to add an additional variable to the backend.
As with all replacement methods, the `length` of `value` has to match
the number of chromatograms represented by the backend. For replacement
of retention time or intensity values we need also to ensure that the
data would be correct after the operation, i.e., that the number of
retention time and intensity values per chromatogram are the identical
and that all retention time and intensity values are numeric. Finally,
we use the
[`validChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
function to ensure that, after replacement, all core chromatogram
variables have the correct data type.

``` r

#' Replace or add a single chromatogram variable.
setReplaceMethod("$", "ChromBackendTest", function(x, name, value) {
    if (length(x) && length(value) != length(x)) {
        stop(
            "length of 'value' needs to match the number of chromatograms ",
            "in object."
        )
    }
    if (name %in% peaksVariables(x)) {
        if (!is.list(value)) {
            stop("The value for peaksData should be a list")
        }
        for (i in seq_along(value)) {
            x@peaksData[[i]][[name]] <- value[[i]]
            validPeaksData(x@peaksData)
        }
    } else {
        x@chromData[, name] <- value
        validChromData(x@chromData)
    }
    x
})
```

We can thus replace an existing chromatogram variable, such as
`msLevel`:

``` r

#' Values before replacement
be$msLevel
```

    ## [1] 1 1

``` r

#' Replace MS levels
be$msLevel <- c(3L, 2L)

#' Values after replacement
be$msLevel
```

    ## [1] 3 2

We can also add a new chromatogram variables:

``` r

#' Add a new chromatogram variable
be$name <- c("A", "B")
be$name
```

    ## [1] "A" "B"

Or also replace intensity values. Below we replace the intensity values
by adding a value of +3 to each.

``` r

#' Replace intensity values
be$msLevel3 <- be$msLevel + 3
be$msLevel3
```

    ## [1] 6 5

#### `peaksData<-`

The `peaksData<-` method should allow to replace the full peaks data
(retention time and intensity value pairs) of all chromatograms in a
backend. As `value`, a `list` of `data.frame` should be provided with
columns names `"rtime"` and `"intensity"`. Because the full peaks data
is provided at once, this method can (and should) support changing also
the number of peaks per chromatogram (while the methods like `rtime<-`
or `$rtime` would not allow).

``` r

#' replacement method for peaks data
setReplaceMethod("peaksData", "ChromBackendTest", function(object, value) {
    if (!is.list(value)) {
        stop("'value' is expected to be a list")
    }
    if (length(object) && length(object) != length(value)) {
        stop("'value' has to be a list with ", length(object), " elements")
    }
    validPeaksData(value)
    object@peaksData <- value
    object
})
```

With this method we can now replace the peaks data of a backend:

``` r

#' Create a list with peaks matrices; our backend has 3 chromatograms
#' thus our `list` has to be of length 3
tmp <- list(
    data.frame(
        rtime = c(12.3, 14.4, 15.4, 16.4),
        intensity = c(200, 312, 354.1, 232)
    ),
    data.frame(
        rtime = c(14.4),
        intensity = c(13.4)
    )
)

be_2 <- be
#' Assign this peaks data to one of our test backends
peaksData(be_2) <- tmp

#' Evaluate that we properly added the peaks data
peaksData(be_2)
```

    ## [[1]]
    ##   rtime intensity
    ## 1  12.3     200.0
    ## 2  14.4     312.0
    ## 3  15.4     354.1
    ## 4  16.4     232.0
    ## 
    ## [[2]]
    ##   rtime intensity
    ## 1  14.4      13.4

### Methods with available default implementations

Default implementations for the `ChromBackend` class are available for a
large number of methods. Thus, any backend extending this class will
automatically inherit these default implementations. Alternative,
class-specific, versions can, but don’t need to be developed. The
default versions are defined in the *R/ChromBackend.R* file, and also
listed in this section. If alternative versions are implemented it
should be ensured that the expected data type is always used for core
chromatogram variables. Use
[`coreChromVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
and
[`corePeaksVariables()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
to list these mandatory data types.

#### `backendParallelFactor()`

The
[`backendParallelFactor()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
function allows a backend to suggest a preferred way it could be split
for parallel processing. The default implementation returns
[`factor()`](https://rdrr.io/r/base/factor.html) (i.e. a `factor` of
length 0) hence not suggesting any specific splitting setup.

``` r

#' Is there a specific way how the object could be best split for
#' parallel processing?
setMethod("backendParallelFactor", "ChromBackend", function(object, ...) {
    factor()
})
```

``` r

backendParallelFactor(be)
```

    ## factor()
    ## Levels:

#### `chromIndex()`

The
[`chromIndex()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
function should return the value for the `"chromIndex"` chromatogram
variable. As a result, an `integer` of length equal to the number of
chromatograms in `object` needs to be returned. The default
implementation is:

``` r

#' get the values for the chromIndex chromatogram variable
setMethod(
    "chromIndex", "ChromBackend",
    function(object, columns = chromVariables(object)) {
        chromData(object, columns = "chromIndex", drop = TRUE)
    }
)
```

The result of calling this method on our test backend:

``` r

chromIndex(be)
```

    ## [1] NA NA

#### `collisionEnergy()`

The
[`collisionEnergy()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
function should return the value for the `"collisionEnergy"`
chromatogram variable. As a result, a `numeric` of length equal to the
number of chromatograms has to be returned. The default implementation
is:

``` r

#' get the values for the collisionEnergy chromatogram variable
setMethod("collisionEnergy", "ChromBackend", function(object) {
    chromData(object, columns = "collisionEnergy", drop = TRUE)
})
```

The result of calling this method on our test backend:

``` r

collisionEnergy(be)
```

    ## [1] NA NA

The default replacement method for the `collisionEnergy` chromatogram
variable is:

``` r

#' Default replacement method for collisionEnergy
setReplaceMethod(
    "collisionEnergy", "ChromBackend", function(object, value) {
        object$collisionEnergy <- value
        object
    }
)
```

This method thus makes use of the `$<-` replacement method we
implemented above. To test this function we replace the collision energy
below.

``` r

#' Replace the collision energy
collisionEnergy(be) <- c(20, 30)
collisionEnergy(be)
```

    ## [1] 20 30

#### `dataOrigin()`, `dataOrigin<-`

The
[`dataOrigin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `dataOrigin<-` methods return or set the value(s) for the
`"dataOrigin"` chromatogram variable. The values for this chromatogram
variable need to be of type `character` (the length equal to the number
of chromatograms). The default implementation for
[`dataOrigin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
is:

``` r

#' Default implementation to access dataOrigin
setMethod("dataOrigin", "ChromBackend", function(object) {
    chromData(object, columns = "dataOrigin", drop = TRUE)
})
```

Below we use this method to access the values of the `dataOrigin`
chromatogram variable.

``` r

#' Access the dataOrigin values
dataOrigin(be)
```

    ## [1] NA NA

The default implementation for `dataOrigin<-` uses, like all defaults
for replacement methods, the `$<-` method:

``` r

#' Default implementation of the `dataOrigin<-` replacement method
setReplaceMethod("dataOrigin", "ChromBackend", function(object, value) {
    object$dataOrigin <- value
    object
})
```

For our backend we can change the values of the `dataOrigin` variable:

``` r

#' Replace the backend's dataOrigin values
dataOrigin(be) <- rep("from somewhere", 2)
dataOrigin(be)
```

    ## [1] "from somewhere" "from somewhere"

#### `intensity()`, `intensity<-`

The
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
and `intensity<-` methods allow to extract or set the intensity values
of the individual chromatograms represented by the backend. The default
for the
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function, which is expected to return a `list` of `numeric` values with
the intensity values of each chromatogram, uses the
[`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
method:

``` r

#' Default method to extract intensity values
setMethod("intensity", "ChromBackend", function(object) {
    if (length(object)) {
        peaksData(object, column = "intensity", drop = TRUE)
    } else {
        list()
    }
})
```

The default replacement method for intensity values uses the `$<-`
method:

``` r

#' Default implementation of the replacement method for intensity values
setReplaceMethod("intensity", "ChromBackend", function(object, value) {
    pd <- peaksData(object)
    if (!is.list(value) || length(pd) != length(value)) {
        stop("'value' should be a list of the same length as 'object'")
    }
    for (i in seq_along(pd)) {
        if (length(value[[i]]) != nrow(pd[[i]])) {
            stop(paste0(
                "Length of 'value[[", i, "]]' does not match ",
                "the number of rows in the intensity of chromatogram: ",
                i, "'"
            ))
        }
    }
    peaksData(object) <- lapply(seq_along(pd), function(i) {
        pd[[i]]$intensity <- value[[i]]
        return(pd[[i]])
    })
    object
})
```

``` r

#' Replace intensity values
intensity(be)[[1]] <- intensity(be)[[1]] + 10
intensity(be)
```

    ## [[1]]
    ## [1]  133.3  163.6 2364.3  253.4
    ## 
    ## [[2]]
    ## [1] 100.0  80.1

#### `isEmpty()`

The [`isEmpty()`](https://rdrr.io/pkg/S4Vectors/man/List-class.html) is
a simple helper function to evaluate whether chromatograms are *empty*,
i.e. have no peaks (retention time and intensity values). It should
return a logical vector of length equal to the number of chromatograms
in the backend with `TRUE` if a chromatogram is empty and `FALSE`
otherwise. The default implementation uses the
[`lengths()`](https://rdrr.io/r/base/lengths.html) method (defined
further below) that returns for each chromatogram the number of
available data points (peaks).

``` r

#' Default implementation for `isEmpty()`
setMethod("isEmpty", "ChromBackend", function(x) {
    lengths(x) == 0L
})
```

``` r

isEmpty(be)
```

    ## [1] FALSE FALSE

#### `isReadOnly()`

As discussed above, backends can also be *read-only*, hence only
allowing to access, but not to change any values (e.g. if the data is
stored in a data base and the connection to this data base does not
support updating or replacing data). In such cases, the default
[`isReadOnly()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
method can be used, which returns always `TRUE`:

``` r

#' Default implementation of `isReadOnly()`
setMethod("isReadOnly", "ChromBackend", function(object) {
    TRUE
})
```

Backends that support changing data values should implement their own
version (like we did above) to return `FALSE` instead:

``` r

isReadOnly(be)
```

    ## [1] FALSE

#### `length()`

The [`length()`](https://rdrr.io/r/base/length.html) method should
return a single `integer` with the total number of chromatograms
available through the backend. The default implementation for this
function is:

``` r

#' Default implementation for `length()`
setMethod("length", "ChromBackend", function(x) {
    nrow(chromData(x, columns = "dataStorage"))
})
```

``` r

length(be)
```

    ## [1] 2

#### `lengths()`

The [`lengths()`](https://rdrr.io/r/base/lengths.html) function should
return the number of data pairs (peaks; retention time or intensity
values) per chromatogram. The result should be an `integer` vector (of
length equal to the number of chromatograms in the backend) with these
counts. The default implementation uses the
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function.

``` r

#' Default implementation for `lengths()`
setMethod("lengths", "ChromBackend", function(x) {
    lengths(intensity(x))
})
```

The number of peaks for our test backend:

``` r

lengths(be)
```

    ## [1] 4 2

#### `msLevel()`, `msLevel<-`

The
[`msLevel()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `msLevel<-` methods should allow extracting and setting the MS level
for the individual chromatograms. MS levels are encoded as `integer`,
thus,
[`msLevel()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
must return an `integer` vector of length equal to the number of
chromatograms of the backend and `msLevel<-` should take/accept such a
vector as input. The default implementations for both methods are shown
below.

``` r

#' Default methods to get or set MS levels
setMethod("msLevel", "ChromBackend", function(object) {
    chromData(object, columns = "msLevel", drop = TRUE)
})
setReplaceMethod("msLevel", "ChromBackend", function(object, value) {
    object$msLevel <- value
    object
})
```

To test these we below replace the MS levels for our test data set and
extract these values again.

``` r

msLevel(be) <- c(1L, 2L)
msLevel(be)
```

    ## [1] 1 2

#### `mz()`, `mz<-`

The
[`mz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `mz<-` methods should allow to extract or set the m/z value for each
chromatogram. The m/z value of a chromatogram is encoded as `numeric`,
thus, the methods are expected to return or accept a `numeric` vector of
length equal to the number of chromatograms. The default implementations
are shown below.

``` r

#' Default implementations to get or set m/z value(s)
setMethod("mz", "ChromBackend", function(object) {
    chromData(object, columns = "mz", drop = TRUE)
})
setReplaceMethod("mz", "ChromBackend", function(object, value) {
    object$mz <- value
    object
})
```

We below set and extract these *target* m/z values.

``` r

mz(be) <- c(314.3, 312.5)
mz(be)
```

    ## [1] 314.3 312.5

#### `mzMax()`, `mzMax<-`

The
[`mzMax()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `mzMax<-` methods should allow to extract or set the upper m/z
boundary for each chromatogram. m/z values are encoded as `numeric`,
thus, the methods are expected to return or accept a `numeric` vector of
length equal to the number of chromatograms. The default implementations
are shown below.

``` r

#' Default implementations to get or set upper m/z limits
setMethod("mzMax", "ChromBackend", function(object) {
    chromData(object, columns = "mzMax", drop = TRUE)
})
setReplaceMethod("mzMax", "ChromBackend", function(object, value) {
    object$mzMax <- value
    object
})
```

Testing these functions by replacing the upper m/z boundary with new
values.

``` r

mzMax(be) <- mz(be) + 0.01
mzMax(be)
```

    ## [1] 314.31 312.51

#### `mzMin(),`mzMin\<-\`

The
[`mzMin()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `mzMin<-` methods should allow to extract or set the lower m/z
boundary for each chromatogram. m/z values are encoded as `numeric`,
thus, the methods are expected to return or accept a `numeric` vector of
length equal to the number of chromatograms. The default implementations
are shown below.

``` r

#' Default methods to get or set the lower m/z boundary
setMethod("mzMin", "ChromBackend", function(object) {
    chromData(object, columns = "mzMin", drop = TRUE)
})

setReplaceMethod("mzMin", "ChromBackend", function(object, value) {
    object$mzMin <- value
    object
})
```

Testing these functions by replacing the lower m/z boundary with new
values.

``` r

mzMin(be) <- mz(be) - 0.01
mzMin(be)
```

    ## [1] 314.29 312.49

#### `precursorMz()`, `precursorMz<-`

The
[`precursorMz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
and `precursorMz<-` methods are expected to get or set the values for
the precursor m/z of each chromatogram (if available). These are encoded
as `numeric` (one value per chromatogram) - and if a value is not
available `NA_real_` should be returned. The default implementations
are:

``` r

#' Default implementations to get or set the precursorMz chrom variable
setMethod("precursorMz", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMz", drop = TRUE)
})
setReplaceMethod("precursorMz", "ChromBackend", function(object, value) {
    object$precursorMz <- value
    object
})
```

Below we set and get the `precursorMz` chromatogram variable for our
backend.

``` r

precursorMz(be) <- c(NA_real_, 123.3)
precursorMz(be)
```

    ## [1]    NA 123.3

#### `precursorMzMax()`, `precursorMzMax<-`

These methods are supposed to allow to get and set the `precursorMzMax`
chromatogram variable. The default implementations are:

``` r

#' Default implementations for `precursorMzMax`
setMethod("precursorMzMax", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMzMax", drop = FALSE)
})
setReplaceMethod("precursorMzMax", "ChromBackend", function(object, value) {
    object$precursorMzMax <- value
    object
})
```

Below we test these functions by setting and extracting the values for
this chromatogram variable.

``` r

precursorMzMax(be) <- precursorMz(be) + 0.1
precursorMzMax(be)
```

    ## [1]    NA 123.4

#### `precursorMzMin()`, `precursorMzMin<-`

These methods are supposed to allow to get and set the `precursorMzMin`
chromatogram variable. The default implementations are:

``` r

#' Default implementations for `precursorMzMin`
setMethod("precursorMzMin", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMzMin", drop = FALSE)
})
setReplaceMethod("precursorMzMin", "ChromBackend", function(object, value) {
    object$precursorMzMin <- value
    object
})
```

Below we test these functions by setting and extracting the values for
this chromatogram variable.

``` r

precursorMzMin(be) <- precursorMz(be) - 0.1
precursorMzMin(be)
```

    ## [1]    NA 123.2

#### `productMz()`, `productMz<-`

These methods are supposed to allow to get and set the `productMz`
chromatogram variable. The default implementations are:

``` r

#' Default implementations for `productMz`
setMethod("productMz", "ChromBackend", function(object) {
    chromData(object, columns = "productMz", drop = TRUE)
})
setReplaceMethod("productMz", "ChromBackend", function(object, value) {
    object$productMz <- value
    object
})
```

Below we test these functions by setting and extracting the values for
this chromatogram variable.

``` r

productMz(be) <- c(123.2, NA_real_)
productMz(be)
```

    ## [1] 123.2    NA

#### `productMzMax()`, `productMzMax<-`

These methods are supposed to allow to get and set the `productMzMax`
chromatogram variable. The default implementations are:

``` r

#' Default implementations for `productMzMax`
setMethod("productMzMax", "ChromBackend", function(object) {
    chromData(object, columns = "productMzMax", drop = FALSE)
})
setReplaceMethod("productMzMax", "ChromBackend", function(object, value) {
    object$productMzMax <- value
    object
})
```

Below we test these functions by setting and extracting the values for
this chromatogram variable.

``` r

productMzMax(be) <- productMz(be) + 0.02
productMzMax(be)
```

    ## [1] 123.22     NA

#### `productMzMin()`, `productMzMin<-`

These methods are supposed to allow to get and set the `productMzMin`
chromatogram variable. The default implementations are:

``` r

#' Default implementations for `productMzMin`
setMethod("productMzMin", "ChromBackend", function(object) {
    chromData(object, columns = "productMzMin", drop = FALSE)
})
setReplaceMethod("productMzMin", "ChromBackend", function(object, value) {
    object$productMzMin <- value
    object
})
```

Below we test these functions by setting and extracting the values for
this chromatogram variable.

``` r

productMzMin(be) <- productMz(be) - 0.2
productMzMin(be)
```

    ## [1] 123  NA

#### `rtime()`, `rtime<-`

The [`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
and `rtime<-` methods allow to get and set the retention times of the
individual chromatograms of the backend. Similar to the method for the
intensity values described above they should return or accept a
`NumericList`, each element being a `numeric` vector with the retention
time values of one chromatogram. The default implementations of these
methods are shown below.

``` r

#' Default methods for `rtime()` and `rtime<-`
setMethod("rtime", "ChromBackend", function(object) {
    if (length(object)) {
        peaksData(object, column = "rtime", drop = TRUE)
    } else {
        list()
    }
})

setReplaceMethod("rtime", "ChromBackend", function(object, value) {
    pd <- peaksData(object)
    if (!is.list(value) || length(pd) != length(value)) {
        stop("'value' should be a list of the same length as 'object'")
    }
    for (i in seq_along(pd)) {
        if (length(value[[i]]) != nrow(pd[[i]])) {
            stop(paste0(
                "Length of 'value[[", i, "]]' does not match ",
                "the number of rows in 'the rtime of chromatogram: ", i, "'"
            ))
        }
    }
    peaksData(object) <- lapply(seq_along(pd), function(i) {
        pd[[i]]$rtime <- value[[i]]
        return(pd[[i]])
    })
    object
})
```

We below test this implementation replacing the retention times of our
example backend by shifting all values by 2 seconds.

``` r

rtime(be)[[1]] <- rtime(be)[[1]] + 2
rtime(be)
```

    ## [[1]]
    ## [1] 14.4 14.8 15.2 16.6
    ## 
    ## [[2]]
    ## [1] 45.1 46.2

#### `split()`

The [`split()`](https://rdrr.io/r/base/split.html) method should split
the backend into a `list` of backends containing subsets of the original
backend. The default implementation uses the default implementation of
[`split()`](https://rdrr.io/r/base/split.html) from R and should work in
most cases. This function uses the `[` method to subset/split the
object.

``` r

#' Default method to split a backend
setMethod("split", "ChromBackend", function(x, f, drop = FALSE, ...) {
    split.default(x, f, drop = drop, ...)
})
```

We below test this by splitting the backend into two subsets.

``` r

split(be, f = c(1, 2, 1))
```

    ## Warning in split.default(x, f, drop = drop, ...): data length is not a multiple
    ## of split variable

    ## $`1`
    ## ChromBackendTest with 1 chromatograms
    ## 
    ## $`2`
    ## ChromBackendTest with 1 chromatograms

## Session information

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] Chromatograms_1.3.3 ProtGenerics_1.45.0 BiocParallel_1.47.0
    ## [4] BiocStyle_2.41.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0         compiler_4.6.1         BiocManager_1.30.27   
    ##  [4] parallel_4.6.1         cluster_2.1.8.3        jquerylib_0.1.4       
    ##  [7] systemfonts_1.3.2      IRanges_2.47.2         textshaping_1.0.5     
    ## [10] yaml_2.3.12            fastmap_1.2.0          R6_2.6.1              
    ## [13] generics_0.1.4         knitr_1.51             BiocGenerics_0.59.10  
    ## [16] htmlwidgets_1.6.4      MASS_7.3-66            bookdown_0.47         
    ## [19] desc_1.4.3             Spectra_1.23.3         bslib_0.12.0          
    ## [22] rlang_1.3.0            cachem_1.1.0           xfun_0.60             
    ## [25] fs_2.1.0               MsCoreUtils_1.25.4     sass_0.4.10           
    ## [28] otel_0.2.0             cli_3.6.6              pkgdown_2.2.1.9000    
    ## [31] digest_0.6.39          MetaboCoreUtils_1.21.1 lifecycle_1.0.5       
    ## [34] clue_0.3-68            S4Vectors_0.51.6       data.table_1.18.4     
    ## [37] evaluate_1.0.5         codetools_0.2-20       ragg_1.5.2            
    ## [40] stats4_4.6.1           rmarkdown_2.31         tools_4.6.1           
    ## [43] htmltools_0.5.9
