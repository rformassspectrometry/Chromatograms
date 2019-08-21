#' @include hidden_aliases.R
NULL

#' @title Chromatographic data backends
#'
#' @aliases ChromBackend-class
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
#' @param mz For `filterMz`: `numeric(2)` defining the lower and upper m/z of
#'     the range to subset `object`. Chromatograms with `mzMin` or `mzMax` within
#'     `mz` are retained.
#'     For `filterPrecursorMz` and `filterProductMz`: `numeric(1)`
#'     with the m/z value to filter the object.
#'
#' @param name For `$` and `$<-`: the name of the chromatogram variable to
#'     return or set.
#'
#' @param object Object extending `ChromBackend`.
#'
#' @param ppm For `filterPrecursorMz` and `filterProductMz`: `numeric(1)`
#'     defining the accepted difference between the provided m/z and the
#'     chromatogrm's (precursor or product) m/z in parts per million.
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
#' - `chromNames`, `chromNames<-`: gets or sets the names for the chromatograms.
#'
#' - `chromVariables`: returns a `character` vector with the available
#'   chromatogram variables available in `object`.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of
#'   chromatograms in `object` with the *data origin* of each. This could e.g.
#'   be the mzML file from which the data was read.
#'
#' - `dataStorage`: gets a `character` of length equal to the number of
#'   chromatograms in `object` with the data storage of each. Note that a
#'   `dataStorage` of `NA_character_` is not supported.
#'
#' - `filterDataOrigin`: filters the object retaining chromatograms matching the
#'   provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   chromatograms to subset.
#'   `filterDataOrigin` should return the data ordered by the provided
#'   `dataOrigin` parameter, i.e. if `dataOrigin = c("2", "1")` was provided,
#'   the chromatograms in the resulting object should be ordered accordingly
#'   (first chromatogram from data origin `"2"` and then from `"1"`).
#'
#' - `filterDataStorage`: filters the object retaining chromatograms matching
#'   the provided `dataStorage`. Parameter `dataStorage` has to be of type
#'   `character` and needs to match exactly the data storage value of the
#'   chromatograms to subset.
#'   `filterDataStorage` should return the data ordered by the provided
#'   `dataStorage` parameter, i.e. if `dataStorage = c("2", "1")` was provided,
#'   the chromatograms in the resulting object should be ordered accordingly
#'   (first chromatogram from data storage `"2"` and then from `"1"`).
#'
#' - `filterMsLevel`: retains chromatograms of MS level `msLevel`.
#'
#' - `filterMz`: retains chromatograms with `mzMin` or `mzMax` within the
#'   provided m/z range.
#'
#' - `filterPrecursorMz`: retains chromatograms with a precursor m/z matching
#'   the provided `mz` accepting also a small difference in m/z which can be
#'   defined by parameter `ppm` (parts per million). With the default
#'   (`ppm = 0`) only chromatograms with m/z identical to `mz` are retained.
#'
#' - `filterProductMz`: retains chromatograms with a product m/z matching the
#'   provided `mz` accepting also a small difference in m/z which can be
#'   defined by parameter `ppm` (parts per million). With the default
#'   (`ppm = 0`) only chromatograms with product m/z identical to `mz` are
#'   retained.
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
#' - `pairs`: returns a `list` of matrices with the data pairs from each
#'   chromatogram. The length of the `list` is equal to the number of
#'   chromatograms in `object`. Each element is a two-column `numeric` `matrix`
#'   with the retention time (first column) and intensity values (second column)
#'   of one chromatogram. For an empty chromatogram a `matrix` with 0 rows and
#'   two columns (named `"rtime"` and `"intensity"`) has to be returned.
#'
#' - `pairs<-`: replaces the chromatogram data of the backend. This function
#'   expects a `list` of two-column matrices in the format returned by `pairs`.
#'   The length of `value` has to be identical to the number of chromatograms
#'   in `object`. Note that just writeable backends support this method.
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
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
#'
#' @importFrom IRanges NumericList
#'
#' @exportClass ChromBackend
NULL

setClass("ChromBackend",
         contains = "VIRTUAL",
         slots = c(
             readonly = "logical",
             version = "character"),
         prototype = prototype(readonly = FALSE, version = "0.1"))

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

#' @exportMethod chromVariables
#'
#' @rdname ChromBackend
setMethod("chromVariables", "ChromBackend", function(object) {
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
setMethod("filterPrecursorMz", "ChromBackend", function(object, mz, ppm) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod filterProductMz
#'
#' @importMethodsFrom ProtGenerics filterProductMz
#'
#' @rdname ChromBackend
setMethod("filterProductMz", "ChromBackend", function(object, mz, ppm, ...) {
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

#' @exportMethod pairs
#'
#' @rdname ChromBackend
setMethod("pairs", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod pairs<-
#'
#' @rdname ChromBackend
setReplaceMethod("pairs", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
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
