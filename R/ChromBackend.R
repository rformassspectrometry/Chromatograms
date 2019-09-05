#' @include hidden_aliases.R
NULL

#' @title Chromatographic data backends
#'
#' @aliases ChromBackend-class ChromBackendDataFrame-class
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
#' The core chromatogram variables are:
#' - `chromIndex`: an `integer` with the index of the chromatogram in the
#'   original source file (e.g. *mzML* file). This is a required variable
#'   for the `ChromBackendMzR` backend but might not be needed or defined for
#'   other backends.
#' - `collisionEnergy`: for SRM data, `numeric` with the collision energy of
#'   the precursor.
#' - `dataOrigin`: optional `character` with the origin of a chromatogram.
#' - `dataStorage`: `character` defining where the data is (currently) stored.
#' - `intensity`: `NumericList` with the intensity values of each chromatogram.
#' - `msLevel`: `integer` defining the MS level of the data.
#' - `mz`: optional `numeric` with the (target) m/z value for the
#'   chromatographic data.
#' - `mzMin`: optional `numeric` with the minimal m/z value in case the data
#'   was extracted from a `Spectra` object.
#' - `mzMax`: optional `numeric` with the maximal m/z value.
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
#' - `rtime`: `NumericList` with the retention times of each chromatogram.
#'
#' @name ChromBackend
#'
#' @param columns For `chromData` accessor: optional `character` with column
#'     names (chromatogram variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param dataOrigin For `filterDataOrigin`: `character` to define which
#'     chromatograms to keep.
#'
#' @param dataStorage For `filterDataStorage`: `character` to define which
#'     chromatograms to keep.
#'
#' @param drop For `[`: ignored.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: ignored.
#'
#' @param msLevel `integer` defining the MS level of the chromatograms to which
#'     the function should be applied. For `filterMsLevel`: the MS level to
#'     which `object` should be subsetted.
#'
#' @param mz For `filterMz`, `filterPrecursorMz` and `filterProductMz`:
#'     `numeric(2)` defining the lower and upper m/z of the range to subset
#'     `object`. Chromatograms with their min m/z (`mzMin`, `precursorMzMin`
#'     and `productMz`, respectively) **or** max m/z within `mz` are retained
#'     (filtering thus returns chromatograms **overlapping** the specified m/z
#'     range).
#'
#' @param name For `$` and `$<-`: the name of the chromatogram variable to
#'     return or set.
#'
#' @param object Object extending `ChromBackend`.
#'
#' @param chromVariables For `selectChromVariables`: `character` with the
#'     names of the chromatogram variables to which the backend should be
#'     subsetted.
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
#' **have** to implement the following methods:
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed
#'
#' - `$`, `$<-`: access or set/add a single chromatogram variable (column) in
#'   the backend.
#'
#' - `as.list`: returns a `list` of matrices with the data pairs from each
#'   chromatogram. The length of the `list` is equal to the number of
#'   chromatograms in `object`. Each element is a two-column `numeric` `matrix`
#'   with the retention time (first column) and intensity values (second column)
#'   of one chromatogram. For an empty chromatogram a `matrix` with 0 rows and
#'   two columns (named `"rtime"` and `"intensity"`) has to be returned.
#'
#' - `backendInitialize`: initialises the backend. This method is
#'   supposed to be called right after creating an instance of the
#'   backend class and should prepare the backend. This method has to ensure
#'   to set the spectra variable `dataStorage` correctly.
#'
#' - `backendMerge`: merges (combines) `ChromBackend` objects into a single
#'   instance. All objects to be merged have to be of the same type.
#'
#' - `chromData`, `chromData<-`: gets or sets general chromatogram
#'   metadata (annotation, also called header).  `chromData` returns
#'   a `DataFrame`, `chromData<-` expects a `DataFrame` with the same number
#'   of rows as there are chromatograms in `object`.
#'
#' - `chromIndex`: returns a `integer` vector with the index of the chromatogram
#'   in the original source file.
#'
#' - `chromNames`, `chromNames<-`: gets or sets the names for the chromatograms.
#'
#' - `chromVariables`: returns a `character` vector with the available
#'   chromatogram variables available in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the collision energy
#'   for the precursor (for SRM data). `collisionEnergy` returns a `numeric` of
#'   length equal to the number of chromatograms in `object`.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of
#'   chromatograms in `object` with the *data origin* of each. This could e.g.
#'   be the mzML file from which the data was read.
#'
#' - `dataStorage`: gets a `character` of length equal to the number of
#'   chromatograms in `object` with the data storage of each. Note that a
#'   `dataStorage` of `NA_character_` is not supported.
#'
#' - `filterDataOrigin`: filters the object retaining chromatograms matching any
#'   of the provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   chromatograms to subset.
#'   `filterDataOrigin` should return the data ordered by the provided
#'   `dataOrigin` parameter, i.e. if `dataOrigin = c("2", "1")` was provided,
#'   the chromatograms in the resulting object should be ordered accordingly
#'   (first chromatogram from data origin `"2"` and then from `"1"`).
#'
#' - `filterDataStorage`: filters the object retaining chromatograms matching
#'   any of the provided `dataStorage`. Parameter `dataStorage` has to be of
#'   type `character` and needs to match exactly the data storage value of the
#'   chromatograms to subset.
#'   `filterDataStorage` should return the data ordered by the provided
#'   `dataStorage` parameter, i.e. if `dataStorage = c("2", "1")` was provided,
#'   the chromatograms in the resulting object should be ordered accordingly
#'   (first chromatogram from data storage `"2"` and then from `"1"`).
#'
#' - `filterMsLevel`: retains chromatograms of MS level `msLevel`.
#'
#' - `filterMz`: retains chromatograms with `mzMin` or `mzMax` within the
#'   provided m/z range (i.e. `mzMax(object) >= mz[1] & mzMin(object) <= mz[2]`
#'   ).
#'
#' - `filterPrecursorMz`: retains chromatograms with a their precursor m/z
#'   window overlapping the provided m/z range `mz` (i.e.
#'   `precursorMzMax(object) >= mz[1] & precursorMzMin(object) <= mz[2]`).
#'
#' - `filterProductMz`: retains chromatograms with a their product m/z
#'   window overlapping the provided m/z range `mz` (i.e.
#'   `productMzMax(object) >= mz[1] & productMzMax(object) <= mz[2]`).
#'
#' - `intensity`: gets the intensity values from the chromatograms. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   chromatogram). The length of the list is equal to the number of
#'   chromatograms in `object`.
#'
#' - `intensity<-`: replaces the intensity values. `value` has to be a `list`
#'   (or [NumericList()]) of length equal to the number of chromatograms and the
#'   number of values within each list element identical to the number of
#'   data pairs in each chromatogram. Note that just writeable backends
#'   support this method.
#'
#' - `isEmpty`: returns a `logical` of length equal to the number of
#'   chromatograms with `TRUE` for chromatograms without any data pairs.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `length`: returns the number of chromatograms in the object.
#'
#' - `lengths`: returns the number of data pairs (rtime-intensity values)
#'   per chromatogram.
#'
#' - `msLevel`: gets the chromatogram's MS level. Returns an `integer`
#'   vector (of length equal to the number of chromatograms) with the MS
#'   level for each chromatogram (or `NA_integer_` if not available).
#'
#' - `mz`,`mz<-`: gets or sets the m/z value of the chromatograms. `mz` returns
#'   a `numeric` of length equal to the number of chromatograms in `object`,
#'   `mz<-` expects a `numeric` of length `length(object)`.
#'
#' - `mzMax`,`mzMax<-`: gets or sets the upper m/z of the mass-to-charge range
#'   from which the chromatogram contains signal (e.g. if the chromatogram
#'   was extracted from MS data in spectra format and a m/z range was provided).
#'   `mzMax` returns a `numeric` of length equal to the number of chromatograms
#'   in `object`, `mzMax` expects a `numeric` of length equal to the number
#'   of chromatograms in `object`.
#'
#' - `mzMin`,`mzMin<-`: gets or sets the lower m/z of the mass-to-charge range
#'   from which the chromatogram contains signal (e.g. if the chromatogram
#'   was extracted from MS data in spectra format and a m/z range was provided).
#'   `mzMin` returns a `numeric` of length equal to the number of chromatograms
#'   in `object`, `mzMin` expects a `numeric` of length equal to the number
#'   of chromatograms in `object`.
#'
#' - `precursorMz`,`precursorMz<-`: gets or sets the (target) m/z of the
#'   precursor (for SRM data). `precursorMz` returns a `numeric` of length
#'   equal to the number of chromatograms in `object`. `precursorMz<-` expects
#'   a `numeric` of length equal to the number of chromatograms.
#'
#' - `productMz`,`productMz<-`: gets or sets the (target) m/z of the
#'   product (for SRM data). `productMz` returns a `numeric` of length
#'   equal to the number of chromatograms in `object`. `productMz<-` expects
#'   a `numeric` of length equal to the number of chromatograms.
#'
#' - `precursorMzMin`,`precursorMzMax`,`productMzMin`, `productMzMax`: get
#'   the lower and upper margin for the precursor or product isolation windows.
#'   These functions might return the value of `productMz` if the respective
#'   minimal or maximal m/z values are not defined in `object`.
#'
#' - `rtime`: gets the retention times from the chromatograms. returns a
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
#' - `selectChromVariables`: reduce `object` retaining only specified
#'   chromatogram variables.
#'
#' @section `ChromBackendDataFrame`, in-memory chromatographic data backend:
#'
#' The `ChromBackendDataFrame` objects keep all chromatographic data in memory.
#' To reduce memory requirement, all chromatogram variables with a single
#' value (e.g. if all are from MS level 1) are internally represented
#' as an [Rle()] object that are converted into the original class (e.g.
#' `integer`) when the column is accessed.
#'
#' New objects can be created with the `ChromBackendDataFrame()`
#' function. The backend can be subsequently initialized with the
#' `backendInitialize` method, taking a `DataFrame` with the chromatographic
#' data as parameter. Suggested columns of this `DataFrame` are:
#'
#' - `"msLevel"`: `integer` with MS levels of the spectra.
#' - `"mz"`: `numeric` with (target) m/z of the spectra.
#' - `"mzMin"`,`"mzMax"`: `numeric` with the lower and upper m/z per
#'   chromatogram.
#' - `"precursorMz"`,`"precursorMzMin"`,`"precursorMzMax"`: m/z value (target,
#'   lower and upper limit) of the precursor ion (for SRM data).
#' - `"productMz"`,`"productMzMin"`,`"productMzMax"`: m/z value (target,
#'   lower and upper limit) of the product ion (for SRM data).
#' - `"dataOrigin"`: `character` defining the *data origin*.
#' - `"dataStorage"`: `character` indicating grouping of spectra in different
#'   e.g. input files. Note that missing values are not supported.
#' - `"precursorMz"`: `numeric` with the m/z value of the precursor.
#' - `"rtime"`: [NumericList()] of `numeric` vectors representing the retention
#'   times for the data values in a chromatogram.
#' - `"intensity"`: [NumericList()] of `numeric` vectors representing the
#'   intensity values for each chromatogram.
#'
#' Additional columns are allowed too.
#'
#' @section `ChromBackendMzR`, on-disc chromatographic data backend:
#'
#' The `ChromBackendMzR` keeps only a limited amount of data in memory,
#' while the chromatographic data (retention time and intensity values) are
#' fetched from the raw files on-demand. This backend uses the `mzR` package
#' for data import and retrieval and hence requires that package to be
#' installed. Also, it can only be used to import and represent data
#' stored in *mzML* files.
#'
#' The `ChromBackendMzR` backend extends the `ChromBackendDataFrame` backend
#' using its `DataFrame` to keep chromatogram variables (except retention time
#' and intensity) in memory.
#'
#' New objects can be created with the `ChromBackendMzR()` function which
#' can be subsequently filled with data by calling `backendInitialize`
#' passing the file names of the input data files with argument `files`.
#'
#' @section Implementation notes:
#'
#' Backends extending `ChromBackend` **must** implement all of its methods
#' (listed above). Developers of new `ChromBackend`s should follow the
#' `ChromBackendDataFrame` implementation.
#'
#' The `ChromBackend` defines the following slots:
#'
#' - `@readonly`: `logical(1)` whether the backend supports writing/replacing
#'   of m/z or intensity values.
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @importFrom IRanges NumericList
#'
#' @exportClass ChromBackend
#'
#' @examples
#'
#' ## In-memory backend
#'
#' ## Create a simple data set with 3 chromatograms
#' df <- data.frame(msLevel = c(1L, 2L, 1L),
#'     mz = c(34.5, 454.2, 123.4))
#' df$rtime <- list(1:4, 2:3, 1:5)
#' df$intensity <- list(c(34, 25, 5, 1), c(5, 2), c(1, 5, 7, 4, 2))
#'
#' be <- backendInitialize(ChromBackendDataFrame(), df)
#' be
#'
#' ## Access variables
#' msLevel(be)
#' be$msLevel
#'
#' rtime(be)
#' be$rtime
#'
#' ## Available variables
#' chromVariables(be)
#'
#' ## Add an additional variable
#' be$new_var <- c("a", "b", "c")
#'
#' chromVariables(be)
#' be$new_var
#'
#' ## Extract the full data
#' chromData(be)
NULL

setClass("ChromBackend",
         contains = "VIRTUAL",
         slots = c(
             readonly = "logical",
             version = "character"),
         prototype = prototype(readonly = FALSE, version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("ChromBackend", function(object) {
    msg <- .valid_chrom_backend_data_storage(dataStorage(object))
    if (length(dataStorage(object)) != length(object))
        msg <- c(msg, "length of object and 'dataStorage' have to match")
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @importFrom methods .valueClassTest validObject
#'
#' @rdname ChromBackend
setMethod("backendInitialize", signature = "ChromBackend",
          definition = function(object, ...) {
              validObject(object)
              object
          })

#' @rdname ChromBackend
setMethod("backendMerge", "list", function(object, ...) {
    backendMerge(object[[1]], object[-1])
})

#' @exportMethod backendMerge
#'
#' @rdname ChromBackend
setMethod("backendMerge", "ChromBackend", function(object, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod chromData
#'
#' @rdname ChromBackend
setMethod("chromData", "ChromBackend",
          function(object, columns = chromVariables(object)) {
              stop("Not implemented for ", class(object), ".")
          })

#' @exportMethod chromData<-
#'
#' @rdname ChromBackend
setReplaceMethod("chromData", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod chromNames
#'
#' @rdname ChromBackend
setMethod("chromNames", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod chromNames<-
#'
#' @rdname ChromBackend
setReplaceMethod("chromNames", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod chromIndex
#'
#' @rdname ChromBackend
setMethod("chromIndex", "ChromBackend",
          function(object, columns = chromVariables(object)) {
              stop("Not implemented for ", class(object), ".")
          })

#' @exportMethod chromVariables
#'
#' @rdname ChromBackend
setMethod("chromVariables", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy
#'
#' @importMethodsFrom ProtGenerics collisionEnergy
#'
#' @rdname ChromBackend
setMethod("collisionEnergy", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod collisionEnergy<-
#'
#' @importMethodsFrom ProtGenerics collisionEnergy<-
#'
#' @rdname ChromBackend
setReplaceMethod("collisionEnergy", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataOrigin
#'
#' @importMethodsFrom ProtGenerics dataOrigin
#'
#' @rdname ChromBackend
setMethod("dataOrigin", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataOrigin<-
#'
#' @importMethodsFrom ProtGenerics dataOrigin<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataOrigin", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname ChromBackend
setMethod("dataStorage", "ChromBackend", function(object) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage<-
#'
#' @importMethodsFrom ProtGenerics dataStorage<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataStorage", "ChromBackend", function(object, value) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
})

#' @exportMethod filterDataOrigin
#'
#' @importMethodsFrom ProtGenerics filterDataOrigin
#'
#' @rdname ChromBackend
setMethod("filterDataOrigin", "ChromBackend", function(object, dataOrigin, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterDataStorage
#'
#' @importMethodsFrom ProtGenerics filterDataStorage
#'
#' @rdname ChromBackend
setMethod("filterDataStorage", "ChromBackend",
          function(object, dataStorage, ...) {
              stop("Not implemented for ", class(object), ".")
          })

#' @exportMethod filterMsLevel
#'
#' @importMethodsFrom ProtGenerics filterMsLevel
#'
#' @rdname ChromBackend
setMethod("filterMsLevel", "ChromBackend", function(object, msLevel) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterMz
#'
#' @importMethodsFrom ProtGenerics filterMz
#'
#' @rdname ChromBackend
setMethod("filterMz", "ChromBackend", function(object, mz, msLevel, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterPrecursorMz
#'
#' @importMethodsFrom ProtGenerics filterPrecursorMz
#'
#' @rdname ChromBackend
setMethod("filterPrecursorMz", "ChromBackend", function(object, mz, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterProductMz
#'
#' @importMethodsFrom ProtGenerics filterProductMz
#'
#' @rdname ChromBackend
setMethod("filterProductMz", "ChromBackend", function(object, mz, ...) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod intensity
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname ChromBackend
setMethod("intensity", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname ChromBackend
setReplaceMethod("intensity", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod isEmpty
#'
#' @rdname ChromBackend
#'
#' @importMethodsFrom S4Vectors isEmpty
setMethod("isEmpty", "ChromBackend", function(x) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod isReadOnly
#'
#' @rdname ChromBackend
setMethod("isReadOnly", "ChromBackend", function(object) {
    object@readonly
})

#' @exportMethod length
#'
#' @rdname ChromBackend
setMethod("length", "ChromBackend", function(x) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod msLevel
#'
#' @importMethodsFrom ProtGenerics msLevel
#'
#' @rdname ChromBackend
setMethod("msLevel", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod msLevel<-
#'
#' @importMethodsFrom ProtGenerics msLevel<-
#'
#' @rdname ChromBackend
setReplaceMethod("msLevel", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mz
#'
#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname ChromBackend
setMethod("mz", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname ChromBackend
setReplaceMethod("mz", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mzMax
#'
#' @rdname ChromBackend
setMethod("mzMax", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("mzMax", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mzMin
#'
#' @rdname ChromBackend
setMethod("mzMin", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod mzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("mzMin", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod as.list
#'
#' @importMethodsFrom BiocGenerics as.list
#'
#' @rdname ChromBackend
setMethod("as.list", "ChromBackend", function(x, ...) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod lengths
#'
#' @rdname ChromBackend
setMethod("lengths", "ChromBackend", function(x) {
    stop("Not implemented for ", class(x), ".")
})

#' @exportMethod precursorMz
#'
#' @importMethodsFrom ProtGenerics precursorMz
#'
#' @rdname ChromBackend
setMethod("precursorMz", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMz<-
#'
#' @importMethodsFrom ProtGenerics precursorMz<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMz", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMzMin
#'
#' @rdname ChromBackend
setMethod("precursorMzMin", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMzMin", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMzMax
#'
#' @rdname ChromBackend
setMethod("precursorMzMax", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod precursorMzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("precursorMzMax", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMz
#'
#' @importMethodsFrom ProtGenerics productMz
#'
#' @rdname ChromBackend
setMethod("productMz", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMz<-
#'
#' @importMethodsFrom ProtGenerics productMz<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMz", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMzMin
#'
#' @rdname ChromBackend
setMethod("productMzMin", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMzMin<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMzMin", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMzMax
#'
#' @rdname ChromBackend
setMethod("productMzMax", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod productMzMax<-
#'
#' @rdname ChromBackend
setReplaceMethod("productMzMax", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod rtime
#'
#' @importMethodsFrom ProtGenerics rtime
#'
#' @rdname ChromBackend
setMethod("rtime", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod rtime<-
#'
#' @importMethodsFrom ProtGenerics rtime<-
#'
#' @rdname ChromBackend
setReplaceMethod("rtime", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod selectChromVariables
#'
#' @rdname ChromBackend
setMethod("selectChromVariables", "ChromBackend",
          function(object, chromVariables = chromVariables(object)) {
              stop("Not implemented for ", class(object), ".")
})

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
