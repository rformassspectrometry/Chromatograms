#' @title Chromatographic MS Data Backends
#'
#' @name ChromBackend
#'
#' @aliases ChromBackend-class
#' @aliases ChromBackendMemory-class
#' @aliases [,ChromBackend-method
#'
#' @description
#'
#' `ChromBackend` is a virtual class that defines what different backends need
#' to provide to be used by the `Chromatograms` package and classes.
#'
#' The backend should provide access to the chromatographic data which mainly
#' consists of (paired) intensity and retention time values. Additional
#' chromatographic metadata such as MS level and precursor and product m/z
#' should also be provided.
#'
#' Through their implementation different backends can be either optimized for
#' minimal memory requirements or performance. Each backend needs to implement
#' data access methods listed in section *Backend functions:* below.
#'
#' And example implementation and more details and descriptions are provided
#' in the *Creating new `ChromBackend` classes for Chromatograms* vignette.
#'
#' Currently available backends are:
#'
#' - `ChromBackendMemory`: This backend stores chromatographic data directly
#'   in memory, making it ideal for small datasets or testing. It can be
#'   initialized with a `data.frame` of chromatographic data via the
#'   `chromData` parameter and a `list` of `data.frame` entries for peaks data
#'   using the `peaksData` parameter. These data can be accessed with the
#'   `chromData()` and `peaksData()` functions.
#'
#'
#' @section Core chromatogram variables:
#'
#' The *core* chromatogram variables are variables (metadata) that can/should
#' be provided by a backend. For each of these variables a value needs to be
#' returned, if none is defined, a missing value (of the correct data type)
#' should be returned. The names of the chromatogram variables in your current
#' chromatogram object are returned with the `chromVariables()` function.
#'
#' For each core chromatogram variable a dedicated access method exists. In
#' contrast to the peaks data described below, a single value should be
#' returned for each chromatogram.
#'
#' The `coreChromVariables()` function returns the core chromatogram variables
#' along with their expected (defined) data type.
#'
#' The core chromatogram variables (in alphabetical order) are:
#'
#' - `chromIndex`: an `integer` with the index of the chromatogram in the
#'   original source file (e.g. *mzML* file).
#' - `collisionEnergy`: for SRM data, `numeric` with the collision energy of
#'   the precursor.
#' - `dataOrigin`: optional `character` with the origin of a chromatogram.
#' - `dataStorage`: `character` defining where the data is (currently) stored.
#' - `msLevel`: `integer` defining the MS level of the data.
#' - `mz`: optional `numeric` with the (target) m/z value for the
#'   chromatographic data.
#' - `mzMin`: optional `numeric` with the lower m/z value of the m/z range in
#'   case the data (e.g. an extracted ion chromatogram EIC) was extracted from
#'   a `Spectra` object.
#' - `mzMax`: optional `numeric` with the upper m/z value of the m/z range.
#' - `precursorMz`: for SRM data, `numeric` with the target m/z of the
#'   precursor (parent).
#' - `precursorMzMin`: for SRM data, optional `numeric` with the lower m/z of
#'   the precursor's isolation window.
#' - `precursorMzMax`: for SRM data, optional `numeric` with the upper m/z of
#'   the precursor's isolation window.
#' - `productMz` for SRM data, `numeric` with the target m/z of the
#'   product ion.
#' - `productMzMin`: for SRM data, optional `numeric` with the lower m/z of
#'   the product's isolation window.
#' - `productMzMax`: for SRM data, optional `numeric` with the upper m/z of
#'   the product's isolation window.
#'
#' @section Core Peaks variables:
#'
#' Similar to the *core* chromatogram variables, *core* peaks variables
#' represent  metadata that should be provided by a backend. Each of these
#' variables should return a value, and if undefined, a missing value (with the
#' appropriate data type) is returned. The number of values for a peaks
#' variable in a single chromatogram can vary, from none to multiple, and may
#' differ between chromatograms.
#'
#' The names of peaks variables in the current chromatogram object can be
#' obtained with the `peaksVariables()` function.
#'
#' Each core peaks variable has a dedicated accessor method.
#'
#'
#' The `corePeaksVariables()` function returns the core peaks variables along
#' with their expected (defined) data type.
#'
#' The core peaks variables, listed in the required order for `peaksData`, are:
#'
#' - `rtime`: A `numeric` vector containing retention time values.
#' - `intensity`: A `numeric` vector containing intensity values.
#'
#' They should be provided for each chromatogram in the backend, **in this order**,
#' No NAs are allowed for the `rtime` values. These characteristics will be
#' checked with the `validPeaksData()` function.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'     for more information.
#'
#' @param columns For `chromData()` accessor: optional `character` with column
#'     names (chromatogram variables) that should be included in the
#'     returned `data.frame`. By default, all columns are returned.
#'
#' @param drop For `chromData()` and `peaksData()`: `logical(1)` default to
#' `FALSE`. If `TRUE`, and one column is called by the user, the method should
#' return a vector (or list of vector for `peaksData()`) of the single column
#' requested.
#'
#' @param f `factor` defining the grouping to split `x`. See [split()].
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[` and `[[`: ignored.
#'
#' @param keep For `filterChromData()`: `logical(1)`
#'        defining whether to keep (`keep = TRUE`) or remove (`keep = FALSE`)
#'        the chromatogram data that match the condition.
#'
#' @param match For `filterChromData()` : `character(1) `
#'        defining whether the condition has to match for all provided
#'        `ranges` (`match = "all"`; the default), or for any of them
#'        (`match = "any"`) for chromatogram data to be retained.
#'
#' @param name For `$` and `$<-`: the name of the chromatogram variable to
#'        return or set.
#'
#' @param object Object extending `ChromBackend`.
#'
#' @param ranges For `filterChromData()` : a `numeric`
#'        vector of paired values (upper and lower boundary) that define the
#'        ranges to filter the `object`. These paired values need to be in the
#'        same order as the `variables` parameter (see below).
#'
#'
#' @param value Replacement value for `<-` methods. See individual
#'        method description or expected data type.
#'
#' @param variables For `filterChromData()`: `character` vector with the names
#'        of the chromatogram variables to filter for. The list of available
#'        chromatogram variables can be obtained with `chromVariables()`.
#'
#' @param x Object extending `ChromBackend`.
#'
#' @param ... Additional arguments.
#'
#' @section Mandatory methods:
#'
#' New backend classes **must** extend the base `ChromBackend` class and
#' implement the following mandatory methods:
#'
#' - `backendInitialize()`: initialises the backend. This method is
#'   supposed to be called right after creating an instance of the
#'   backend class and should prepare the backend.
#'   Parameters can be defined freely for each backend, depending on what is
#'   needed to initialize the backend.
#'   This method has to ensure to set the chromtogram variable `dataStorage`
#'   correctly.
#'
#' - `backendBpparam()`: returns the parallel processing setup supported by
#'   the backend class. This function can be used by any higher
#'   level function to evaluate whether the provided parallel processing
#'   setup (or the default one returned by `bpparam()`) is supported
#'   by the backend. Backends not supporting parallel processing (e.g.
#'   because they contain a connection to a database that can not be
#'   shared across processes) should extend this method to return only
#'   `SerialParam()` and hence disable parallel processing for (most)
#'   methods and functions. See also `backendParallelFactor()` for a
#'   function to provide a preferred splitting of the backend for parallel
#'   processing.
#'
#' - `backendParallelFactor()`: returns a `factor` defining an optimal
#'   (preferred) way how the backend can be split for parallel processing
#'   used for all peak data accessor or data manipulation functions.
#'   The default implementation returns a factor of length 0 (`factor()`)
#'   providing thus no default splitting. `backendParallelFactor()` for
#'   `ChromBackendMzR` on the other hand returns `factor(dataStorage(object))`
#'   hence suggesting to split the object by data file.
#'
#' - `chromData()`, `chromData<-`: gets or sets general chromatogram metadata
#'   (annotation). `chromData()` returns a `data.frame`, `chromData<-` expects
#'   a `data.frame` with the same number of rows as there are chromatograms in
#'   `object`. Read-only backends might not need to implement the
#'   replacement method `chromData<-` (unless some internal caching mechanism
#'   could be used). `chromData()` should be implemented with the parameter
#'   `drop` set to `FALSE` as default. With `drop = FALSE` the method should
#'   return a `data.frame` even if only one column is called. If `drop = TRUE`
#'   is specified, the output will be a vector of the single column requested.
#'   New backends should be implemented such as if empty, the method returns a
#'   `data.frame` with 0 rows and the columns defined by `chromVariables()`.
#'   By default, the function *should* return at minimum the coreChromVariables,
#'   even if NAs.
#'
#' - `peaksData()`: returns a `list` of `data.frame` with the data
#'   (e.g. retention time - intensity pairs) from each chromatogram. The length
#'   of the `list` is equal to the number of chromatograms in `object`. For an
#'   empty chromatogram a `data.frame` with 0 rows and two columns (named
#'   `"rtime"` and `"intensity"`) has to be returned. The optional parameter
#'   `columns`, if supported by the backend allows to define which peak
#'   variables should be returned in each array. As default (minimum) columns
#'   `"rtime"` and `"intensity"` have to be provided. `peaksData()` should be
#'   implemented with the parameter `drop` set to `FALSE` as default.  With
#'   `drop = FALSE` the method should return a `data.frame` even if only one
#'   column is called. If `drop = TRUE`  is specified, the output will be a
#'   vector of the single column requested.
#'
#' - `peaksData<-` replaces the peak data (retention time and intensity values)
#'   of the backend. This method expects a `list` of two-dimensional arrays
#'   (`data.frame`) with columns representing the peak variables.
#'   All existing peaks data are expected to be replaced with these new values.
#'   The length of the `list` has to match the number of chromatogram of `object`.
#'   Note that only writeable backends need to support this method.
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed. This method should be implemented as to support empty integer.
#'
#' - `$`, `$<-`: access or set/add a single chromatogram variable (column) in
#'   the backend.
#'
#' - `backendMerge()`: merges (combines) `ChromBackend` objects into a single
#'   instance. All objects to be merged have to be of the same type.
#'
#' @section Optional methods with default implementations:
#'
#' Additional methods that might be implemented, but for which default
#' implementations are already present are:
#'
#' - `[[`
#'
#' - `backendParallelFactor()`: returns a `factor` defining an optimal
#'   (preferred) way how the backend can be split for parallel processing
#'   used for all *peak* data accessor or data manipulation functions.
#'   The default implementation returns a factor of length 0 (`factor()`)
#'   providing thus no default splitting.
#'
#' - `chromIndex()`: returns an `integer` vector with the index of the
#'   chromatograms in the original source file.
#'
#' - `chromVariables()`: returns a `character` vector with the
#'   available chromatogram variables (columns, fields or attributes)
#'   available in `object`. Variables listed by this function are expected to
#'   be returned (if requested) by the `chromData()` function.
#'
#' - `collisionEnergy()`, `collisionEnergy<-`: gets or sets the collision energy
#'   for the precursor (for SRM data). `collisionEnergy()` returns a `numeric`
#'   of length equal to the number of chromatograms in `object`.
#'
#' - `dataOrigin()`, `dataOrigin<-`: gets or sets the *data origin* variable.
#'   `dataOrigin()` returns a `character` of length equal to the number of
#'   chromatograms, `dataOrigin<-` expects a `character` of length equal
#'   `length(object)`.
#'
#' - `dataStorage()`, `dataStorage<-`: gets or sets the *data storage* variable.
#'   `dataStorage()` returns a `character` of length equal to the number of
#'   chromatograms in `object`, `dataStorage<- ` expects a `character` of
#'   length equal `length(object)`. Note that missing values (`NA_character_`)
#'   are not supported for `dataStorage()`.
#'
#' - `filterChromData()`: filters any numerical chromatographic data variables
#'   based on the provided numerical `ranges`. The method should return a
#'   `ChromBackend` object with the chromatograms that match the condition. This
#'   function will results in an object with less chromatogram than the original.
#'
#' - `intensity()`: gets the intensity values from the chromatograms. Returns
#'   a `list` of `numeric` vectors (intensity values for each
#'   chromatogram). The length of the list is equal to the number of
#'   chromatograms in `object`.
#'
#' - `intensity<-`: replaces the intensity values. `value` has to be a `list`
#'   of length equal to the number of chromatograms and the number of values
#'   within each list element identical to the number of data pairs in each
#'   chromatogram. Note that just writeable backends need to support this
#'   method.
#'
#' - `isReadOnly()`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data. Defaults to FALSE.
#'
#' - `isEmpty()`: returns a `logical` of length equal to the number of
#'   chromatograms with `TRUE` for chromatograms without any data pairs.
#'
#' - `length()`: returns the number of chromatograms in the object.
#'
#' - `lengths()`: returns the number of data pairs (retention time and intensity
#'   values) per chromatogram.
#'
#' - `msLevel()`: gets the chromatogram's MS level. Returns an `integer`
#'   vector (of length equal to the number of chromatograms) with the MS
#'   level for each chromatogram (or `NA_integer_` if not available).
#'
#' - `mz()`,`mz<-`: gets or sets the m/z value of the chromatograms. `mz()`
#'   returns a `numeric` of length equal to the number of chromatograms in `
#'   object`, `mz<-` expects a `numeric` of length `length(object)`.
#'
#' - `mzMax()`,`mzMax<-`: gets or sets the upper m/z of the mass-to-charge
#'   range from which a chromatogram contains signal (e.g. if the chromatogram
#'   was extracted from MS data in spectra format and a m/z range was provided).
#'   `mzMax()` returns a `numeric` of length equal to the number of
#'   chromatograms in `object`, `mzMax<-` expects a `numeric` of length equal
#'   to the number of chromatograms in `object`.
#'
#' - `mzMin()`,`mzMin<-`: gets or sets the lower m/z of the mass-to-charge range
#'   from which a chromatogram contains signal (e.g. if the chromatogram
#'   was extracted from MS data in spectra format and a m/z range was provided).
#'   `mzMin()` returns a `numeric` of length equal to the number of
#'   chromatograms in `object`, `mzMin<-` expects a `numeric` of length equal
#'   to the number of chromatograms in `object`.
#'
#' - `peaksVariables()`: lists the available data variables for the
#'   chromatograms. Default peak variables are `"rtime"` and `"intensity"`
#'   (which all backends need to support and provide), but some backends
#'   might provide additional variables.
#'   Variables listed by this function are expected to be returned (if
#'   requested) by the `peaksData()` function.
#'
#' - `precursorMz()`,`precursorMz<-`: gets or sets the (target) m/z of the
#'   precursor (for SRM data). `precursorMz()` returns a `numeric` of length
#'   equal to the number of chromatograms in `object`. `precursorMz<-` expects
#'   a `numeric` of length equal to the number of chromatograms.
#'
#' - `precursorMzMin()`,`precursorMzMax()`,`productMzMin()`, `productMzMax()`:
#'   gets the lower and upper margin for the precursor or product isolation
#'   windows. These functions might return the value of `productMz()` if the
#'   respective minimal or maximal m/z values are not defined in `object`.
#'
#' - `productMz()`,`productMz<-`: gets or sets the (target) m/z of the
#'   product (for SRM data). `productMz()` returns a `numeric` of length
#'   equal to the number of chromatograms in `object`. `productMz<-` expects
#'   a `numeric` of length equal to the number of chromatograms.
#'
#' - `rtime()`: gets the retention times from the chromatograms. returns a
#'   [NumericList()] of `numeric` vectors (retention times for each
#'   chromatogram). The length of the returned list is equal to the number of
#'   chromatograms in `object`.
#'
#' - `rtime<-`: replaces the retention times. `value` has to be a `list` (or
#'   [NumericList()]) of length equal to the number of chromatograms and the
#'   number of values within each list element identical to the number of
#'   data pairs in each chromatogram. Note that just writeable backends support
#'   this method.
#'
#' - `split()`: splits the backend into a `list` of backends (depending on
#'   parameter `f`). The default method for `ChromBackend` uses
#'   [split.default()], thus backends extending `ChromBackend` don't
#'   necessarily need to implement this method.
#'
#' @section Implementation notes:
#'
#' Backends extending `ChromBackend` **must** implement all of its methods
#' (listed above). A guide to create new backend classes is provided as a
#' dedicated vignette. Additional information and an example for a backend
#' implementation is provided in the respective vignette.
#'
#' @author Johannes Rainer, Philippine Louail
#'
#' @md
#'
#' @importFrom IRanges NumericList
#'
#' @exportClass ChromBackend
#'
#' @examples
#'
#' ## Create a simple backend implementation
#' ChromBackendDummy <- setClass("ChromBackendDummy",
#'     contains = "ChromBackend")
NULL

setClass("ChromBackend",
         contains = "VIRTUAL",
         slots = c(
             version = "character"),
         prototype = prototype(readonly = FALSE, version = "0.1"))

#' @importMethodsFrom S4Vectors [ [<-
#' @exportMethod [
#'
#' @rdname ChromBackend
setMethod("[", "ChromBackend", function(x, i, j, ..., drop = FALSE) {
    stop("Not implemented for ", class(x), ".")
})

#' @importMethodsFrom S4Vectors $ $<-
#'
#' @exportMethod $
#'
#' @rdname ChromBackend
setMethod("$", "ChromBackend", function(x, name) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod $<-
#'
#' @rdname ChromBackend
setReplaceMethod("$", "ChromBackend", function(x, name, value) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod backendMerge
#'
#' @importMethodsFrom ProtGenerics backendMerge
#'
#' @rdname ChromBackend
setMethod("backendMerge", "ChromBackend", function(object, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod chromData
#'
#' @rdname ChromBackend
setMethod("chromData", "ChromBackend",
          function(object, columns = chromVariables(object), drop = FALSE) {
              stop("Not implemented for ", class(object), ".")
          })

#' @exportMethod chromData<-
#'
#' @rdname ChromBackend
setReplaceMethod("chromData", "ChromBackend",
                 function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod peaksData
#'
#' @importMethodsFrom ProtGenerics peaksData
#'
#' @rdname ChromBackend
setMethod("peaksData", "ChromBackend",
    function(object, columns = c("rtime", "intensity"), drop = FALSE) {
        stop("Not implemented for ", class(object), ".")
    })

#' @exportMethod peaksData<-
#'
#' @importMethodsFrom ProtGenerics peaksData<-
#'
#' @rdname ChromBackend
setReplaceMethod("peaksData", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

################################################################################
## Methods with default implementations below.

#' @rdname ChromBackend
#'
#' @importMethodsFrom S4Vectors [[ [[<-
#'
#' @export
setMethod("[[", "ChromBackend", function(x, i, j, ...) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the chromatogram ",
             "variable to access.")
    if (!missing(j))
        stop("'j' is not supported.")
    do.call("$", list(x, i))
})

#' @rdname ChromBackend
#'
#' @export
setReplaceMethod("[[", "ChromBackend", function(x, i, j, ..., value) {
    if (!is.character(i))
        stop("'i' is supposed to be a character defining the chromatogram ",
             "variable to replace or create.")
    if (!missing(j))
        stop("'j' is not supported.")
    do.call("$<-", list(x, i, value))
})

#' @importMethodsFrom ProtGenerics backendBpparam
#'
#' @exportMethod backendBpparam
#'
#' @rdname ChromBackend
#'
#' @export
setMethod("backendBpparam", signature = "ChromBackend",
          function(object, BPPARAM = bpparam()) {
              BPPARAM
          })

#' @exportMethod backendInitialize
#'
#' @importFrom methods .valueClassTest validObject
#'
#' @importMethodsFrom ProtGenerics backendInitialize
#'
#' @rdname ChromBackend
setMethod("backendInitialize", signature = "ChromBackend",
          definition = function(object, ...) {
              validObject(object)
              object
          })

#' @rdname ChromBackend
#'
#' @importMethodsFrom ProtGenerics backendParallelFactor
#'
#' @exportMethod backendParallelFactor
setMethod("backendParallelFactor", "ChromBackend", function(object, ...) {
    factor()
})

#' @rdname ChromBackend
#'
#' @importMethodsFrom ProtGenerics backendMerge
setMethod("backendMerge", "list", function(object, ...) {
    backendMerge(object[[1]], object[-1])
})

#' @exportMethod chromIndex
#'
#' @rdname ChromBackend
setMethod("chromIndex", "ChromBackend",
          function(object) {
              chromData(object, columns = "chromIndex", drop = TRUE)
              })

#' @exportMethod chromIndex<-
#'
#' @rdname ChromBackend
setReplaceMethod("chromIndex", "ChromBackend", function(object, value) {
    object$chromIndex <- value
    object
})

#' @exportMethod chromVariables
#'
#' @rdname ChromBackend
setMethod("chromVariables", "ChromBackend", function(object) {
    names(coreChromVariables())
})

#' @exportMethod collisionEnergy
#'
#' @importMethodsFrom ProtGenerics collisionEnergy
#'
#' @rdname ChromBackend
setMethod("collisionEnergy", "ChromBackend", function(object) {
    chromData(object, columns = "collisionEnergy", drop = TRUE)
})

#' @exportMethod collisionEnergy<-
#'
#' @importMethodsFrom ProtGenerics collisionEnergy<-
#'
#' @rdname ChromBackend
setReplaceMethod("collisionEnergy", "ChromBackend", function(object, value) {
    object$collisionEnergy <- value
    object
})

#' @exportMethod dataOrigin
#'
#' @importMethodsFrom ProtGenerics dataOrigin
#'
#' @rdname ChromBackend
setMethod("dataOrigin", "ChromBackend", function(object) {
    chromData(object, columns = "dataOrigin", drop = TRUE)
})

#' @exportMethod dataOrigin<-
#'
#' @importMethodsFrom ProtGenerics dataOrigin<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataOrigin", "ChromBackend", function(object, value) {
    object$dataOrigin <- value
    object
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname ChromBackend
setMethod("dataStorage", "ChromBackend", function(object) {
    chromData(object, columns = "dataStorage", drop = TRUE)
})

#' @exportMethod dataStorage<-
#'
#' @importMethodsFrom ProtGenerics dataStorage<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataStorage", "ChromBackend", function(object, value) {
    object$dataStorage <- value
    object
})

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname ChromBackend
setMethod("intensity", "ChromBackend", function(object) {
    if (length(object)) {
        peaksData(object, column = "intensity", drop = TRUE)
    } else list()
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname ChromBackend
setReplaceMethod("intensity", "ChromBackend", function(object, value) {
    pd <- peaksData(object)
    if (!is.list(value) || length(pd) != length(value))
        stop("'value' should be a list of the same length as 'object'")
    for (i in seq_along(pd)) {
        if (length(value[[i]]) != nrow(pd[[i]]))
            stop(paste0("Length of 'value[[", i, "]]' does not match ",
                       "the number of rows in the intensity of chromatogram: ",
                       i, "'"))
    }
    peaksData(object) <- lapply(seq_along(pd), function(i) {
        pd[[i]]$intensity <- value[[i]]
        return(pd[[i]])
    })
    object
})

#' @exportMethod isEmpty
#'
#' @rdname ChromBackend
#'
#' @importMethodsFrom S4Vectors isEmpty
setMethod("isEmpty", "ChromBackend", function(x) {
    lengths(x) == 0L
})

#' @exportMethod isReadOnly
#'
#' @importMethodsFrom ProtGenerics isReadOnly
#'
#' @rdname ChromBackend
setMethod("isReadOnly", "ChromBackend", function(object) FALSE)

#' @exportMethod length
#'
#' @rdname ChromBackend
setMethod("length", "ChromBackend", function(x) {
    nrow(chromData(x, columns = "dataStorage"))
})

#' @exportMethod lengths
#'
#' @rdname ChromBackend
setMethod("lengths", "ChromBackend", function(x) {
    lengths(intensity(x))
})

#' @exportMethod msLevel
#'
#' @importMethodsFrom ProtGenerics msLevel
#'
#' @rdname ChromBackend
setMethod("msLevel", "ChromBackend", function(object) {
    chromData(object, columns = "msLevel", drop = TRUE)
})

#' @exportMethod msLevel<-
#'
#' @importMethodsFrom ProtGenerics msLevel<-
#'
#' @rdname ChromBackend
setReplaceMethod("msLevel", "ChromBackend", function(object, value) {
    object$msLevel <- value
    object
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname ChromBackend
setMethod("mz", "ChromBackend", function(object) {
    chromData(object, columns = "mz", drop = TRUE)
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname ChromBackend
setReplaceMethod("mz", "ChromBackend", function(object, value) {
    object$mz <- value
    object
})

#' @exportMethod mzMax
#'
#' @rdname ChromBackend
setMethod("mzMax", "ChromBackend", function(object) {
    chromData(object, columns = "mzMax", drop = TRUE)
})

#' @exportMethod mzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("mzMax", "ChromBackend", function(object, value) {
    object$mzMax <- value
    object
})

#' @exportMethod mzMin
#'
#' @rdname ChromBackend
setMethod("mzMin", "ChromBackend", function(object) {
    chromData(object, columns = "mzMin", drop = TRUE)
})

#' @exportMethod mzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("mzMin", "ChromBackend", function(object, value) {
    object$mzMin <- value
    object
})

#' @exportMethod peaksData
#'
#' @importMethodsFrom ProtGenerics peaksVariables
#'
#' @rdname ChromBackend
setMethod("peaksVariables", "ChromBackend", function(object) {
    names(corePeaksVariables())
})

#' @exportMethod precursorMz
#'
#' @importMethodsFrom ProtGenerics precursorMz
#'
#' @rdname ChromBackend
setMethod("precursorMz", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMz", drop = TRUE)
    })

#' @exportMethod precursorMz<-
#'
#' @importMethodsFrom ProtGenerics precursorMz<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMz", "ChromBackend", function(object, value) {
    object$precursorMz <- value
    object
})

#' @exportMethod precursorMzMax
#'
#' @rdname ChromBackend
setMethod("precursorMzMax", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMzMax", drop = TRUE)
})

#' @exportMethod precursorMzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMzMax", "ChromBackend", function(object, value) {
    object$precursorMzMax <- value
    object
})

#' @exportMethod precursorMzMin
#'
#' @rdname ChromBackend
setMethod("precursorMzMin", "ChromBackend", function(object) {
    chromData(object, columns = "precursorMzMin", drop = TRUE)
})

#' @exportMethod precursorMzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMzMin", "ChromBackend", function(object, value) {
    object$precursorMzMin <- value
    object
})

#' @exportMethod productMz
#'
#' @importMethodsFrom ProtGenerics productMz
#'
#' @rdname ChromBackend
setMethod("productMz", "ChromBackend", function(object) {
    chromData(object, columns = "productMz", drop = TRUE)
})

#' @exportMethod productMz<-
#'
#' @importMethodsFrom ProtGenerics productMz<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMz", "ChromBackend", function(object, value) {
    object$productMz <- value
    object
})

#' @exportMethod productMzMax
#'
#' @rdname ChromBackend
setMethod("productMzMax", "ChromBackend", function(object) {
    chromData(object, columns = "productMzMax", drop = TRUE)
})

#' @exportMethod productMzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMzMax", "ChromBackend", function(object, value) {
    object$productMzMax <- value
    object
})

#' @exportMethod productMzMin
#'
#' @rdname ChromBackend
setMethod("productMzMin", "ChromBackend", function(object) {
    chromData(object, columns = "productMzMin", drop = TRUE)
})

#' @exportMethod productMzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMzMin", "ChromBackend", function(object, value) {
    object$productMzMin <- value
    object
})

#' @exportMethod reset
#'
#' @rdname ChromBackend
setMethod("reset", "ChromBackend", function(object) {
    object
})

#' @exportMethod rtime
#'
#' @importMethodsFrom ProtGenerics rtime
#'
#' @rdname ChromBackend
setMethod("rtime", "ChromBackend", function(object) {
    if (length(object)) {
        peaksData(object, column = "rtime", drop = TRUE)
    } else list()
})

#' @exportMethod rtime<-
#'
#' @importMethodsFrom ProtGenerics rtime<-
#'
#' @rdname ChromBackend
setReplaceMethod("rtime", "ChromBackend", function(object, value) {
    pd <- peaksData(object)
    if (!is.list(value) || length(pd) != length(value))
        stop("'value' should be a list of the same length as 'object'")
    for (i in seq_along(pd)) {
        if (length(value[[i]]) != nrow(pd[[i]])) {
            stop(paste0("Length of 'value[[", i, "]]' does not match ",
            "the number of rows in 'the rtime of chromatogram: ", i, "'"))
        }
    }
    peaksData(object) <- lapply(seq_along(pd), function(i) {
        pd[[i]]$rtime <- value[[i]]
        return(pd[[i]])
    })
    object
})

#' @exportMethod split
#'
#' @importMethodsFrom S4Vectors split
#'
#' @rdname ChromBackend
setMethod("split", "ChromBackend", function(x, f, drop = FALSE, ...) {
    split.default(x, f, drop = drop, ...)
})

################################################################################
## Filter functions

#' @exportMethod filterChromData
#'
#' @rdname ChromBackend
#'
#' @importFrom MsCoreUtils between
#'
setMethod("filterChromData", "ChromBackend",
          function(object, variables = character(),
                   ranges = numeric(), match = c("any", "all"),
                   keep = TRUE) {
              if (!length(variables) || !length(ranges))
                  return(object)
              if (!is.numeric(ranges))
                  stop("filterChromData only support filtering for numerical ",
                       "'variables'")
              match <- match.arg(match)
              if (is.character(variables)){
                  if(!all(variables %in% chromVariables(object)))
                      stop("One or more values passed with parameter ",
                           "'variables' are not available as chromatogram ",
                           "variables in object. Use the 'chromVariables()' ",
                           "function to list possible values.")
              } else
                  stop("The 'variables' parameter needs to be of type ",
                       "'character'.")
              if (length(variables) != length(ranges) / 2)
                  stop("Length of 'ranges' needs to be twice the length of ",
                       "the parameter 'variables' and define the lower ",
                       "and upper bound for values of each chromatogram ",
                       "variable defined with parameter 'variables'.")
              query <- chromData(object, columns = variables)
              idx <- .filter_ranges(query, ranges, match)
              if (keep) return(object[idx])
              else {
                  if (length(idx)) return(object[-idx])
                  else return(object)
              }
          })

#' @exportMethod filterPeaksData
#'
#' @rdname ChromBackend
#'
#' @description
#' Filter the peak data based on the provided ranges for the given variables.
#'
#' @note This function replaces the peaksData() of the input object. Therefore
#' backend with `readOnly == TRUE` (i.e. ChromBackendmzR) will need to have a
#' carefully implemented `peaksData(object) <-` method.
#' First thought are to implement a way  to "cache" the results in a slot of
#' the object (maybe an actual `@peaksData` slot or  a `cacheData` slot).
#'
#' The important things to be aware of are:
#'  - This slot should only be used temporarily as to not overload the memory.
#'  - Maybe inspiring on the `MsBackendCached` class and method.
#'  - E.g. not storing the data but some sort of indices ?
#'  - Storing the data in a  temporary file and reading it back when needed ?
#'
#' @export
setMethod("filterPeaksData", "ChromBackend",
          function(object, variables = character(),
                   ranges = numeric(), match = c("any", "all"),
                   keep = TRUE) {
              if (!length(ranges) || !length(variables))
                  return(object)
              if (!is.numeric(ranges))
                    stop("filterPeaksData only support filtering for numerical ",
                       "peak variables")
              match <- match.arg(match)
              if (is.character(variables)) {
                  if (!all(variables %in% peaksVariables(object)))
                      stop("One or more values passed with parameter ",
                           "'variables' are not available as peaks variables in ",
                           "object. Use the 'peaksVariables()' function to list ",
                           "possible values.")
              } else
                  stop("The 'variables' parameter needs to be of type ",
                       "'character'.")
              if (length(variables) != length(ranges) / 2)
                  stop("Length of 'ranges' needs to be twice the length of the ",
                       "parameter 'variables' and define the lower and upper ",
                       "bound for values of each peak variable defined with ",
                       "parameter 'variables'.")
              if (keep) sel_fun <- function(z, idx) z[idx, , drop = FALSE]
              else sel_fun <- function(z, idx) {
                  if (idx == 0) return(z) # check this
                  else return(z[-idx, , drop = FALSE]) }
              peaksData(object) <- lapply(peaksData(object), function(pd) {
                  sel_fun(pd, .filter_ranges(pd[, variables, drop = FALSE],
                                             ranges, match))
              })
              object
          })


