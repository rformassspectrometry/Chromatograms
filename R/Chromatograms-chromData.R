#' @include Chromatograms.R hidden_aliases.R

#' @title Chromatographic Peaks Metadata.
#'
#' @name chromData
#'
#' @aliases chromData
#' @aliases msLevel
#' @aliases chromIndex
#' @aliases collisionEnergy
#' @aliases dataOrigin
#' @aliases mz
#' @aliases mzMin
#' @aliases mzMax
#' @aliases precursorMz
#' @aliases precursorMzMin
#' @aliases precursorMzMax
#' @aliases productMz
#' @aliases productMzMin
#' @aliases productMzMax
#' @aliases filterChromData
#' @aliases chromVariables
#' @aliases chromData<-
#' @aliases msLevel<-
#' @aliases chromIndex<-
#' @aliases collisionEnergy<-
#' @aliases dataOrigin<-
#' @aliases mz<-
#' @aliases mzMin<-
#' @aliases mzMax<-
#' @aliases precursorMz<-
#' @aliases precursorMzMin<-
#' @aliases precursorMzMax<-
#' @aliases productMz<-
#' @aliases productMzMin<-
#' @aliases productMzMax<-
#' @aliases chromVariables<-
#'
#' @description
#'
#' As explained in the [`Chromatograms`] class documentation, the
#' `Chromatograms` object is a container for chromatogram data that includes
#' chromatographic peaks data (*retention time* and related intensity values,
#' also referred to as *peaks data variables* in the context of
#' `Chromatograms`) and metadata of individual chromatograms (so called
#'  *chromatograms variables*).
#'
#' The *chromatograms variables* information can be accessed using the
#' `chromData()` function. it is also possible to access specific
#' chromatograms variables using `$`.
#'
#' `chromData` can be accessed, replaced but also filtered/subsetted. Refer to
#' the sections below for more details.
#'
#' @param columns A `character` vector of chromatograms variables to extract.
#'
#' @param drop A `logical` indicating whether to drop dimensions when
#'        extracting a single variable.
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
#' @param ranges For `filterChromData()` : a `numeric`
#'        vector of paired values (upper and lower boundary) that define the
#'        ranges to filter the `object`. These paired values need to be in the
#'        same order as the `variables` parameter (see below).
#'
#' @param object A [Chromatograms] object.
#'
#' @param value replacement value for `<-` methods. See individual
#'        method description or expected data type.
#'
#' @param variables For `filterChromData()`: `character` vector with the names
#'        of the chromatogram variables to filter for. The list of available
#'        chromatogram variables can be obtained with `chromVariables()`.
#'
#' @param x A [Chromatograms] object.
#'
#' @section Chromatograms variables and accessor functions:
#'
#' The following chromatograms variables are guaranteed to be provided by a
#' `Chromatograms` object and to be accessible with either the `chromData()` or
#'  a specific function named after the variables names:
#'
#' - `chromIndex`: an `integer` with the index of the chromatogram in the
#'   original source file (e.g. *mzML* file).
#' - `collisionEnergy`: for SRM data, `numeric` with the collision energy of
#'   the precursor.
#' - `dataOrigin`: optional `character` with the origin of the data.
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
#' @section Filter Chromatograms variables:
#'
#' Functions that filter `Chromatograms` based on chromatograms variables
#' (i.e, `chromData` ) will remove chromatographic data that do not meet the
#' specified conditions. This means that if a chromatogram is filtered out, its
#' corresponding `chromData` and `peaksData` will be removed from the object
#' immediately.
#'
#' The available functions to filter chromatogram data are:
#'
#' - `filterChromData()`: Filters numerical chromatographic data variables
#'   based on the provided numerical `ranges`. The method returns a
#'   `Chromatograms` object containing only the chromatograms that match the
#'   specified conditions. This function results in an object with fewer
#'   chromatograms than the original.
#'
#'
#' @seealso [Chromatograms] for a general description of the `Chromatograms`
#'          object.
#'          [peaksData] for a general description of the chromatographic peaks
#'          data available in the object, as well as how to access, replace and
#'          subset them.
#'          [processingQueue] for more information on the queuing
#'          of processings and parallelization for larger dataset processing.
#' @md
#'
#' @author Philippine Louail
#'
#' @examples
#'
#' # Create a Chromatograms object
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     chromIndex = c(1L, 2L, 3L)
#' )
#'
#' be <- backendInitialize(new("ChromBackendMemory"), chromData = cdata)
#'
#' chr <- Chromatograms(be)
#'
#' # Access chromatograms variables
#' chromData(chr)
#'
#' # Access specific chromatograms variables
#' chromData(chr, columns = "msLevel")
#'
#' msLevel(chr)
#'
#' # Replace chromatograms variables
#' msLevel(chr) <- c(1L, 2L, 2L)
#'
#' # Filter chromatograms variables
#' filterChromData(chr,
#'     variables = "msLevel", ranges = c(1L, 1L),
#'     keep = FALSE
#' )
#'
NULL

#' @rdname chromData
setMethod(
    "chromData", "Chromatograms",
    function(object, columns = chromVariables(object), drop = FALSE) {
        chromData(object@backend, columns = columns, drop = drop)
    }
)

#' @rdname chromData
setReplaceMethod("chromData", "Chromatograms", function(object, value) {
    chromData(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("chromVariables", "Chromatograms", function(object) {
    chromVariables(object@backend)
})

#' @rdname chromData
setMethod("chromIndex", "Chromatograms", function(object) {
    chromIndex(object@backend)
})

#' @rdname chromData
setReplaceMethod("chromIndex", "Chromatograms", function(object, value) {
    chromIndex(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("collisionEnergy", "Chromatograms", function(object) {
    collisionEnergy(object@backend)
})

#' @rdname chromData
setReplaceMethod("collisionEnergy", "Chromatograms", function(object, value) {
    collisionEnergy(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("dataOrigin", "Chromatograms", function(object) {
    dataOrigin(object@backend)
})

#' @rdname chromData
setReplaceMethod("dataOrigin", "Chromatograms", function(object, value) {
    dataOrigin(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("msLevel", "Chromatograms", function(object) msLevel(object@backend))

#' @rdname chromData
setReplaceMethod("msLevel", "Chromatograms", function(object, value) {
    msLevel(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("mz", "Chromatograms", function(object) mz(object@backend))

#' @rdname chromData
setReplaceMethod("mz", "Chromatograms", function(object, value) {
    mz(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("mzMax", "Chromatograms", function(object) mzMax(object@backend))

#' @rdname chromData
setReplaceMethod("mzMax", "Chromatograms", function(object, value) {
    mzMax(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("mzMin", "Chromatograms", function(object) mzMin(object@backend))

#' @rdname chromData
setReplaceMethod("mzMin", "Chromatograms", function(object, value) {
    mzMin(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("length", "Chromatograms", function(x) length(x@backend))

#' @rdname chromData
setMethod("precursorMz", "Chromatograms", function(object) {
    precursorMz(object@backend)
})

#' @rdname chromData
setReplaceMethod("precursorMz", "Chromatograms", function(object, value) {
    precursorMz(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("precursorMzMin", "Chromatograms", function(object) {
    precursorMzMin(object@backend)
})

#' @rdname chromData
setReplaceMethod("precursorMzMin", "Chromatograms", function(object, value) {
    precursorMzMin(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("precursorMzMax", "Chromatograms", function(object) {
    precursorMzMax(object@backend)
})

#' @rdname chromData
setReplaceMethod("precursorMzMax", "Chromatograms", function(object, value) {
    precursorMzMax(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("productMz", "Chromatograms", function(object) {
    productMz(object@backend)
})

#' @rdname chromData
setReplaceMethod("productMz", "Chromatograms", function(object, value) {
    productMz(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("productMzMin", "Chromatograms", function(object) {
    productMzMin(object@backend)
})

#' @rdname chromData
setReplaceMethod("productMzMin", "Chromatograms", function(object, value) {
    productMzMin(object@backend) <- value
    object
})

#' @rdname chromData
setMethod("productMzMax", "Chromatograms", function(object) {
    productMzMax(object@backend)
})

#' @rdname chromData
setReplaceMethod("productMzMax", "Chromatograms", function(object, value) {
    productMzMax(object@backend) <- value
    object
})

#' @rdname chromData
setMethod(
    "filterChromData", "Chromatograms",
    function(object,
    variables = character(), ranges = numeric(),
    match = c("any", "all"), keep = TRUE) {
        object@backend <- filterChromData(object@backend,
            variables = variables,
            ranges = ranges,
            match = match,
            keep = keep
        )
        object
    }
)
