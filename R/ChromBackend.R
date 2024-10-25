#' @title Chromatographic MS Data Backends
#'
#' @name ChromBackend
#'
#' @aliases ChromBackend-class
#' @aliases ChromBackendMemory-class
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
#' - `ChromBackendMemory`: This backend is used to store chromatographic data in
#'   memory. It is particularly useful for small datasets or for testing
#'   purposes. The backend can be initialized by inputing a `data.frame`  with
#'   the chromatographic data in the `chromData` parameter and a `list` of
#'   `data.frame` with the peaks data in the `peaksData` parameter. Both of
#'   these information can be retrieved using the accessor functions
#'   `chromData()` and `peaksData()`.
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
#' For each core chromatogram variable a dedicated access method exists.
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
#' In a similar fashion as the core chromatogram variables, the core peaks
#' variables are variables (metadata) that can/should be provided by a backend.
#' For each of these variables a value needs to be returned, if none is defined,
#' a missing value (of the correct data type) should be returned.
#' The names of the peaks variables in your current chromatogram object are
#' returned with the `peaksVariables()` function.
#'
#' For each core peaks variable a dedicated access method exists.
#'
#' The `corePeaksVariables()` function returns the core peaks variables along
#' with their expected (defined) data type.
#'
#' The core peaks variables (in alphabetical order) are:
#'
#' - `rtime`: `numeric` with the retention time values.
#' - `intensity`: `numeric` with the intensity values.
#'
#' They should be provided for each chromatogram in the backend, **in this order**,
#' No NAs are allowed for the `rtime` values. These characteristics will be
#' checked with the `validPeaksData()` function.
#'
#' @param chromData For `backendInitialize()` of a `ChromBackendMemory` backend,
#'     a `data.frame` with the chromatographic data. If not provided
#'     (or if empty), a default `data.frame` with the core chromatographic
#'     variables will be created.
#'
#' @param columns For `chromData()` accessor: optional `character` with column
#'     names (chromatogram variables) that should be included in the
#'     returned `data.frame`. By default, all columns are returned.
#'
#' @param dataOrigin For `filterDataOrigin()`: `character` to define which
#'     chromatograms to keep.
#'
#' @param dataStorage For `filterDataStorage()`: `character` to define which
#'     chromatograms to keep.
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
#' @param j For `[`: ignored.
#'
#' @param msLevel `integer` defining the MS level of the chromatograms to which
#'     the function should be applied. For `filterMsLevel()`: the MS level to
#'     which `object` should be subsetted.
#'
#' @param mz For `filterMzValues()`: `numeric` with the m/z values of
#'     chromatograms to keep. All chromatograms with their `mz` chromatogram
#'     variable matching any of the values provided with this parameter are
#'     retained. Parameters `ppm` and `tolerance` allow relaxed matching.
#'     For `filterMzRange()`: `numeric(2)` defining the lower and upper
#'     boundary of the m/z range. Chromatograms with their `mz` chromatogram
#'     variable within this range are retained.
#'
#' @param name For `$` and `$<-`: the name of the chromatogram variable to
#'     return or set.
#'
#' @param object Object extending `ChromBackend`.
#'
#' @param peaksData For `backendInitialize()` of a `ChromBackendMemory` backend,
#'     a `list` of `data.frame` with the peaks data. If not provided (or if
#'     empty), a default `list` of empty `data.frame` with the core peaks
#'     variables will be created. The length of the list should match the number
#'     of chromatograms in the `chromData` parameter.
#'
#' @param ppm For `filterMzValues()`: m/z-relative acceptable difference (in
#'     parts-per-million) for m/z values to be considered *matching*.
#'
#' @param tolerance For `filterMzValues()`: largest acceptable absolute
#'     difference in m/z values to consider them *matching*.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `ChromBackend`.
#'
#' @param ... Additional arguments.
#'
#' @section Backend functions:
#'
#' New backend classes **must** extend the base `ChromBackend` class and
#' implement the following mandatory methods:
#'
#' **Phili: I have now implemented a lot of default methods, so should I maybe
#'  remove these from the list below ? or keep them and say that if not
#'  implemented with the new backend the developper has to ensure it satifies
#'  the required described below when ran with the default method. Example:
#'  most accessor method, chromIndex, collisionEnergy, dataOrigin, dataStorage,
#'   msLevel, mz BUT also the filter functions, filterDataOrigin, ...**
#'
#' - `backendInitialize()`: initialises the backend. This method is
#'   supposed to be called right after creating an instance of the
#'   backend class and should prepare the backend.
#'   Parameters can be defined freely for each backend, depending on what is
#'   needed to initialize the backend.
#'   This method has to ensure to set the spectra variable `dataStorage`
#'   correctly.
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
#' - `chromVariables()`: returns a `character` vector with the
#'   available chromatogram variables (columns, fields or attributes)
#'   available in `object`. Vairables listed by this function are expected to
#'   be returned (if requested) by the `chromData()` function.
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
#'   The length of the `list` has to match the number of spectra of `object`.
#'   Note that only writeable backends need to support this method.
#'
#' - `peaksVariables()`: lists the available data variables for the
#'   chromatograms. Default peak variables are `"rtime"` and `"intensity"`
#'   (which all backends need to support and provide), but some backends
#'   might provide additional variables.
#'   Variables listed by this function are expected to be returned (if
#'   requested) by the `peaksData()` function.
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed.
#'
#' - `$`, `$<-`: access or set/add a single chromatogram variable (column) in
#'   the backend.
#'
#' - `backendMerge()`: merges (combines) `ChromBackend` objects into a single
#'   instance. All objects to be merged have to be of the same type.
#'
#' Additional methods that might be implemented, but for which default
#' implementations are already present are:
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
#' - `collisionEnergy()`, `collisionEnergy<-`: gets or sets the collision energy
#'   for the precursor (for SRM data). `collisionEnergy()` returns a `numeric`
#'   of length equal to the number of chromatograms in `object`.
#'
#' - `dataOrigin()`, `dataOrigin<-`: gets or sets the *data origin* variable.
#'   `dataOrigin()` returns a `character` of length equal to the number of
#'   chromatograms, `dataOrigin<-` expects a `character` of length equal
#'   `length(object)`. Note that missing values (`NA_character_`)
#'   are not supported for `dataOrigin()`.
#'
#' - `dataStorage()`, `dataStorage<-`: gets or sets the *data storage* variable.
#'   `dataStorage()` returns a `character` of length equal to the number of
#'   chromatograms in `object`, `dataStorage<- ` expects a `character` of
#'   length equal `length(object)`. Note that missing values (`NA_character_`)
#'   are not supported for `dataStorage()`.
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
#'   only* or does allow also to write/update data.
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
#'
#' Filter methods:
#'
#' - `filterDataOrigin()`: filters the object retaining chromatograms matching
#'   any of the provided `dataOrigin`. Parameter `dataOrigin` has to be of
#'   type `character` and needs to match exactly the data origin value of the
#'   chromatograms to subset.
#'   `filterDataOrigin()` should return the data ordered by the provided
#'   `dataOrigin` parameter, i.e. if `dataOrigin = c("2", "1")` was provided,
#'   the chromatograms in the resulting object should be ordered accordingly
#'   (first chromatogram from data origin `"2"` and then from `"1"`).
#'
#' - `filterDataStorage()`: filters the object retaining chromatograms matching
#'   any of the provided `dataStorage`. Parameter `dataStorage` has to be
#'   of type `character` and needs to match exactly the data storage value of
#'   the chromatograms to subset.
#'   `filterDataStorage()` should return the data ordered by the provided
#'   `dataStorage` parameter, i.e. if `dataStorage = c("2", "1")` was
#'   provided, the chromatograms in the resulting object should be ordered
#'   accordingly (first chromatogram from data storage `"2"` and then from
#'   `"1"`).
#'
#' - `filterMsLevel()`: retains chromatograms of MS level `msLevel()`.
#'
#' - `filterMzRange()`: retains chromatograms with their m/z within the
#'   provided m/z range.
#'
#' - `filterMzValues()`: retains chromatograms with their m/z matching any of
#'   the provided m/z values (given the provided acceptable differences defined
#'   by parameters `tolerance` and `ppm`.
#'
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

#' @exportMethod [
#'
#' @rdname ChromBackend
setMethod("[", "ChromBackend", function(x, i, j, ..., drop = FALSE) {
    stop("Not implemented for ", class(x), ".")
})

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

#' @exportMethod chromVariables
#'
#' @rdname ChromBackend
setMethod("chromVariables", "ChromBackend", function(object) {
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

#' @exportMethod peaksData
#'
#' @importMethodsFrom ProtGenerics peaksVariables
#'
#' @rdname ChromBackend
setMethod("peaksVariables", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

################################################################################
## Methods with default implementations below.

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
    chromData(object)$chromIndex <- value
    object
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
    chromData(object)$collisionEnergy <- value
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
    chromData(object)$dataOrigin <- value
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
    chromData(object)$dataStorage <- value
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
        if (length(value[[i]]) != nrow(pd[[i]])) {
            stop(paste0("Length of 'value[[", i, "]]' does not match ",
                       "the number of rows in the intensity of chromatogram: ",
                       i, "'"))
        }
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
setMethod("isReadOnly", "ChromBackend", function(object) {
    TRUE
})

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
    chromData(object)$msLevel <- value
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
    chromData(object)$mz <- value
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
    chromData(object)$mzMax <- value
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
    chromData(object)$mzMin <- value
    object
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
    chromData(object)$precursorMz <- value
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
    chromData(object)$precursorMzMax <- value
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
    chromData(object)$precursorMzMin <- value
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
    chromData(object)$productMz <- value
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
    chromData(object)$productMzMax <- value
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
    chromData(object)$productMzMin <- value
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
## Filter functions TODO ADD MORE!

#' @exportMethod filterDataOrigin
#'
#' @importMethodsFrom ProtGenerics filterDataOrigin
#'
#' @rdname ChromBackend
setMethod("filterDataOrigin", "ChromBackend",
          function(object, dataOrigin = character(), ...) {
              if (length(dataOrigin)) {
                  object <- object[dataOrigin(object) %in% dataOrigin]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object), dataOrigin))]
                  else object
              } else object
          })

#' @exportMethod filterDataStorage
#'
#' @importMethodsFrom ProtGenerics filterDataStorage
#'
#' @rdname ChromBackend
setMethod("filterDataStorage", "ChromBackend",
          function(object, dataStorage = character()) {
              if (length(dataStorage)) {
                  object <- object[dataStorage(object) %in% dataStorage]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object), dataStorage))]
                  else object
              } else object
          })

#' @exportMethod filterMsLevel
#'
#' @importMethodsFrom ProtGenerics filterMsLevel
#'
#' @rdname ChromBackend
setMethod("filterMsLevel", "ChromBackend",
          function(object, msLevel = integer()) {
              if (length(msLevel)) {
                  object[msLevel(object) %in% msLevel]
              } else object
          })

#' @exportMethod filterMzRange
#'
#' @importMethodsFrom ProtGenerics filterMzRange
#'
#' @importFrom MsCoreUtils between
#'
#' @rdname ChromBackend
setMethod("filterMzRange", "ChromBackend", function(object, mz = numeric(),
                                                    ...) {
    if (length(mz)) {
        mz <- range(mz)
        keep <- which(between(mz(object), mz))
        object[keep]
    } else object
})

#' @exportMethod filterMzValues
#'
#' @importMethodsFrom ProtGenerics filterMzValues
#'
#' @rdname ChromBackend
setMethod("filterMzValues", "ChromBackend",
          function(object, mz = numeric(), ppm = 20, tolerance = 0, ...) {
              if (length(mz)) {
                  object[.values_match_mz(precursorMz(object), mz = mz,
                                          ppm = ppm, tolerance = tolerance)]
              } else object
          })

## TODO:
## - generic filterRanges, filterValues methods as in Spectra.
## - filterRtime hsould be an important one no ?
