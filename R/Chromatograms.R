#' @include ChromBackend.R ChromBackendMemory.R
NULL

#' @title The Chromatograms class to manage and access chromatographic data
#'
#' @name Chromatograms
#'
#' @aliases Chromatograms-class Chromatograms
#' @aliases [,Chromatograms-method
#' @aliases [<-,Chromatograms-method
#' @aliases [[,Chromatograms-method
#' @aliases [[<-,Chromatograms-method
#'
#' @description
#' The `Chromatograms` class encapsules chromatographic data and related
#' metadata. The chromatographic data is represented by a *backend* extending
#' the virtual [ChromBackend] class which provides the raw data to the
#' `Chromatograms` object. Different backends and their properties are
#' decribed in the [ChromBackend] class documentation.
#'
#' @section Creation of objects:
#'
#' `Chromatograms` objects can be created using the `Chromatograms()`
#' construction function. Either by providing a `ChromBackend` object or by
#' providing a `Spectra` object. The `Spectra` object will be used to generate
#' a `Chromatograms` object with a backend of class [`ChromBackendSpectra`].
#'
#' @section Data stored in a `Chromatograms` object:
#'
#' The `Chromatograms` object is a container for chromatographic data, which
#' includes peaks data (*retention time* and related intensity values, also
#' referred to as *peaks data variables* in the context of `Chromatograms`) and
#' metadata of individual chromatogram (so called *chromatograms variables*).
#' While a core set of chromatograms variables (the
#' `coreChromatogramsVariables()`) and peaks data variables (the
#' `corePeaksVariables()`) are guaranteed to be provided by a `Chromatograms`,
#' it is possible to add arbitrary variables to a `Chromatograms` object.
#'
#' The `Chromatograms` object is designed to contain chromatographic data of a
#' (large) set of chromatograms. The data is organized *linearly* and can be
#' thought of a list of chromatograms, i.e. each element in the `Chromatograms`
#' is one chromatogram.
#'
#' The *chromatograms variables* information in the `Chromatograms` object can
#' be accessed using the `chromData()` function. Specific chromatograms
#' variables can be accessed by either precising the `"columns"` parameter in
#' `chromData()` or using `$`. `chromData` can be accessed, replaced but
#' also filtered/subsetted. Refer to the [chromData] documentation for more
#' details.
#'
#' The *peaks data variables* information in the `Chromatograms` object can be
#' accessed using the `peaksData()` function. Specific peaks variables can be
#' accessed by either precising the `"columns"` parameter in `peaksData()` or
#' using `$`. `peaksData` can be accessed, replaced but also
#' filtered/subsetted. Refer to the [peaksData] documentation for more details.
#'
#' @section Processing of `Chromatograms` objects:
#'
#' Functions that process the chromatograms data in some ways can be applied to
#' the object either directly or by using the `processingQueue` mechanism. The
#' `processingQueue` is a list of processing steps that are stored within the
#' object and only applied when needed. This was created so that the data can be
#' processed in a single step and is very useful for larger datasets. This is
#' even more true as this processing queue will call function that can be
#' applied on the data in a chunk-wise manner. This allows for parallel
#' processing of the data and reduces the memory demand. To read more about the
#' `processingQueue`, and how to parallelize your processes, see the
#' [processingQueue] documentation.
#'
#' @param BPPARAM Parallel setup configuration. See [BiocParallel::bpparam()]
#'        for more information.
#'
#' @param backend [ChromBackend] object providing the raw data for the
#'        `Chromatograms` object.
#'
#' @param chromData For `Chromatograms()` build from a `Spectra` object backend,
#'        a `data.frame` with the chromatographic data. If not provided
#'        (or if empty), a default `data.frame` with the core chromatographic
#'        variables will be created.
#'
#' @param drop For `[`: `logical(1)` default to `FALSE`.
#'
#' @param f `factor` defining the grouping to split the `Chromatograms` object.
#'
#' @param factorize.by A `character` vector with the names of the variables in
#'        the `Spectra` object and the `chromData` slot that should be used
#'        to factorize the `Spectra` object data to generate the
#'        chromatographic data.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[` and `[[`: ignored.
#'
#' @param name A `character` string specifying the name of the variable to
#'        access.
#'
#' @param object A [Chromatograms] object.
#'
#' @param processingQueue [list] a list of processing steps (i.e. functions) to
#'        be applied to the chromatographic data. The processing steps are
#'        applied in the order they are listed in the `processingQueue`.
#'
#' @param summarize.method For Chromatograms created with a `Spectra` object:
#'        A `character` vector with the name of the function to be used to
#'        summaries the spectra data intensity. The available methods are "sum"
#'        and "max". The default is "sum".
#'
#' @param value The value to replace the variable with.
#'
#' @param x A [Chromatograms] object.
#'
#' @param ... Additional arguments.
#'
#' @md
#'
#' @exportClass Chromatograms
#'
#' @return Refer to the individual function description for information on
#' the return value.
#'
#' @seealso [chromData] for a general description of the chromatographic
#'          metadata available in the object, as well as how to access, replace
#'          and subset them.
#'          [peaksData] for a general description of the chromatographic peaks
#'          data available in the object, as well as how to access, replace and
#'          subset them.
#'          [processingQueue] for more information on the queuing of
#'          processings and parallelization for larger dataset.
#'
#' @examples
#'
#' library(MsBackendMetaboLights)
#' library(Spectra)
#'
#' ## Create a Chromatograms object from a Spectra object.
#' be <- backendInitialize(MsBackendMetaboLights(),
#'     mtblsId = "MTBLS39",
#'     filePattern = c("63B.cdf")
#' )
#' s <- Spectra(be)
#' be <- backendInitialize(new("ChromBackendSpectra"), s)
#' chr <- Chromatograms(be)
#'
#' ## Subset
#' chr[1:2]
#'
#' ## access a specific variables
#' chr[["msLevel"]]
#' chr$msLevel
#'
#' ## Replace data of a specific variable
#' chr$msLevel <- c(2L, 2L, 2L)
#'
#' ## Can re factorize the data
#' chr <- factorize(chr)
#'
#' ## Can also change the backend into memory
#' chr <- setBackend(chr, ChromBackendMemory())
#'
#' chr
#'
NULL

setClassUnion("ChromBackendOrMissing", c("ChromBackend", "missing"))

#' The Chromatograms class
#'
#' The `Chromatograms` class is a container for chromatographic data.
#'
#' @slot backend [ChromBackend] the *backend* object providing the raw data for
#'       the `Chromatograms` object.
#'
#' @slot processingQueue [list] a list of processing steps to be applied to the
#'       `Chromatograms`. Each element in the list is a function that
#'       processes the data. The processing steps are applied in
#'       the order they are listed in the `processingQueue`.
#'
#' @slot processing [character] a character vector with the names of the
#'       processing steps that have been applied to the `Chromatograms` object.
#'       This is mainly used in the "show" method to display the processing
#'       steps that have been applied to the `Chromatograms` object.
#'
#' @slot processingChunkSize [numeric(1)] the number of chromatograms to be
#'       processed in a single chunk. This is useful for processing large
#'       data sets in smaller chunks to avoid memory issues.
#'
#' @slot version [character(1)] the version of the `Chromatograms` object.
#'
#' @noRd
setClass("Chromatograms",
         slots = c(
             backend = "ChromBackend",
             processingQueue = "list",
             processing = "character",
             processingChunkSize = "numeric",
             version = "character"
         ),
         prototype = prototype(
             version = "0.1",
             processingChunkSize = Inf,
             processingQueue = list(),
             processing = character()
         )
)

setValidity("Chromatograms", function(object) {
    msg <- character()
    if (!is(.backend(object), "ChromBackend")) {
        msg <- ("backend must be a ChromBackend object")
    }
    if (!is.numeric(processingChunkSize(object)) ||
        length(processingChunkSize(object)) != 1) {
        msg <- c(msg, "processingChunkSize must be a numeric value")
    }
    msg <- c(msg, .valid_processing_queue(.processingQueue(object)))
    if (length(msg)) {
        msg
    } else {
        TRUE
    }
})

#' @rdname Chromatograms
#' @export
setMethod(
    "Chromatograms", "ChromBackendOrMissing",
    function(object = ChromBackendMemory(),
             processingQueue = list(), ...) {
        if (missing(object)) {
            object <- ChromBackendMemory()
        }
        new("Chromatograms",
            backend = object,
            processingQueue = processingQueue, ...
        )
    }
)

#' @rdname Chromatograms
#' @importFrom methods new
#' @export
setMethod(
    "Chromatograms", "Spectra",
    function(object, summarize.method = c("sum", "max"),
             chromData = data.frame(),
             factorize.by = c("msLevel", "dataOrigin"), ...) {
        bd <- backendInitialize(ChromBackendSpectra(),
                                spectra = object,
                                factorize.by = factorize.by,
                                chromData = chromData,
                                summarize.method = summarize.method
        )
        new("Chromatograms",
            backend = bd,
            processingQueue = list(), ...
        )
    }
)


#' @rdname hidden_aliases
#'
#' @param object A [Chromatograms] object.
#' @importMethodsFrom methods show
#' @importFrom utils capture.output
#'
#' @exportMethod show
setMethod(
    "show", "Chromatograms",
    function(object) {
        cat("Chromatographic data (", class(object)[1L], ") with ",
            length(.backend(object)), " chromatograms in a ",
            class(.backend(object)), " backend:\n",
            sep = ""
        )
        if (length(.backend(object))) {
            txt <- capture.output(show(.backend(object)))
            cat(txt[-1], sep = "\n")
        }
        if (length(.processingQueue(object))) {
            cat(
                "Lazy evaluation queue:", length(.processingQueue(object)),
                "processing step(s)\n"
            )
        }
        lp <- length(.processing(object))
        if (lp) {
            lps <- .processing(object)
            if (lp > 3) {
                lps <- lps[seq_len(3)]
                lps <- c(lps, paste0(
                    "...", lp - 3,
                    " more processings. ",
                    "Use 'processingLog' to list all."
                ))
            }
            cat("Processing:\n", paste(lps, collapse = "\n "), "\n")
        }
    }
)

#' @rdname Chromatograms
#'
#' @note
#' This needs to be discussed, if we want for example to be able to set a
#' a backend to `ChromBackendMzR` we need to implement backendInitialize()
#' better. = Support peaksData and chromData as arguments AND have a way to
#' write .mzml files (which we do not have for chromatographic data).
#'
#' @importMethodsFrom ProtGenerics setBackend
#'
#' @exportMethod setBackend
setMethod(
    "setBackend", c("Chromatograms", "ChromBackend"),
    function(object, backend, f = processingChunkFactor(object),
             BPPARAM = SerialParam(), ...) {
        backend_class <- class(.backend(object))
        BPPARAM <- backendBpparam(.backend(object), BPPARAM)
        BPPARAM <- backendBpparam(backend, BPPARAM)
        if (!supportsSetBackend(backend)) {
            stop(class(backend), " does not support 'setBackend'")
        }
        if (!length(f) || length(levels(f)) == 1 || !length(object)) {
            bd_new <- backendInitialize(backend,
                                        peaksData = peaksData(object),
                                        chromData = chromData(object)
            )
        } else {
            bd_new <- bplapply(
                split(.backend(object), f = f),
                function(z, ...) {
                    backendInitialize(backend,
                                      peaksData = peaksData(z),
                                      chromData = chromData(z),
                                      BPPARAM = SerialParam()
                    )
                }, ...,
                BPPARAM = BPPARAM
            )
            bd_new <- backendMerge(bd_new)
        }
        if (colnames(chromData(bd_new)) %in% c("rtMin", "rtMax"))
            chromData(bd_new) <- chromData(bd_new)[, - which(colnames(chromData(bd_new)) %in% c("rtMin", "rtMax"))]
        object@backend <- bd_new
        object@processing <- .logging(
            object@processing,
            "Switch backend from ",
            backend_class, " to ",
            class(.backend(object))
        )
        object
    }
)

#' @rdname Chromatograms
#' @export
setMethod("$", signature = "Chromatograms", function(x, name) {
    .backend(x)[[name]]
})

#' @rdname Chromatograms
#' @export
setReplaceMethod("$", signature = "Chromatograms", function(x, name, value) {
    x@backend[[name]] <- value
    x
})

#' @rdname Chromatograms
#' @importFrom methods slot<-
#' @importFrom MsCoreUtils i2index
#' @export
setMethod("[", "Chromatograms", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j)) {
        stop("Subsetting 'Chromatograms' by columns is not (yet) supported")
    }
    if (missing(i)) {
        return(x)
    }
    slot(x, "backend", check = FALSE) <- extractByIndex(
        .backend(x), i2index(i, length(x))
    )
    x
})

#' @rdname Chromatograms
#' @export
setMethod("[[", "Chromatograms", function(x, i, j, ...) {
    if (!is.character(i)) {
        stop(
            "'i' is supposed to be a character defining the chromatogram or ",
            "peak variable to access."
        )
    }
    if (!missing(j)) {
        stop("'j' is not supported.")
    }
    if (!(i %in% peaksVariables(x)) && !(i %in% chromVariables(x))) {
        stop("No variable '", i, "' available")
    } else {
        do.call("[[", list(.backend(x), i))
    }
})

#' @rdname Chromatograms
#'
#' @export
setReplaceMethod("[[", "Chromatograms", function(x, i, j, ..., value) {
    if (!is.character(i)) {
        stop(
            "'i' is supposed to be a character defining the chromatogram ",
            "or peak variable to replace or create."
        )
    }
    if (!(i %in% peaksVariables(x)) && !(i %in% chromVariables(x))) {
        stop("No variable '", i, "' available")
    }
    if (!missing(j)) {
        stop("'j' is not supported.")
    }
    x@backend <- do.call("[[<-", list(.backend(x), i = i, value = value))
    x
})

#' @rdname Chromatograms
#' @export
setMethod(
    "factorize", "Chromatograms",
    function(object, factorize.by = c("msLevel", "dataOrigin"), ...) {
        object@backend <- factorize(.backend(object), ...)
        object
    }
)
