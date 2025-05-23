#' @include helpers.R
#' @include hidden_aliases.R
#' @include ChromBackend.R
NULL

#' @title Improved in-memory Chromatographic data backend
#'
#' @name ChromBackendMemory
#'
#'
#' @description
#' `ChromBackendMemory`: This backend stores chromatographic data directly
#' in memory, making it ideal for small datasets or testing. It can be
#' initialized with a `data.frame` of chromatographic data via the `chromData`
#' parameter and a `list` of `data.frame` entries for peaks data using the
#' `peaksData` parameter. These data can be accessed with the `chromData()` and
#' `peaksData()` functions.
#'
#' @param chromData For `backendInitialize()` of a `ChromBackendMemory`
#'        backend, a `data.frame` with the chromatographic data. If not
#'        provided (or if empty), a default `data.frame` with the core
#'        chromatographicvariables will be created.
#'
#' @param object  A `ChromBackendMemory` object.
#'
#' @param peaksData For `backendInitialize()` of a `ChromBackendMemory`
#'        backend, a `list` of `data.frame` with the peaks data. If not
#'        provided (or if empty), a default `list` of empty `data.frame` with
#'        the core peaks variables will be created. The length of the list
#'        should match the number of chromatograms in the `chromData`
#'        parameter.
#'
#' @param ... Additional parameters to be passed.
#'
#' @author Philippine Louail
#'
#' @exportClass ChromBackendMemory
#'
#' @examples
#'
#' ## Create a ChromBackendMemory object
#' cbm <- ChromBackendMemory()
#'
#' ## Initialize the ChromBackendMemory object with a data.frame of
#' ## chromatographic data  and a list of data.frame of peaks data
#' cdata <- data.frame(
#'     msLevel = c(1L, 1L, 1L),
#'     mz = c(112.2, 123.3, 134.4),
#'     chromIndex = c(1L, 2L, 3L)
#' )
#'
#' pdata <- list(
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     ),
#'     data.frame(
#'         rtime = c(45.1, 46.2),
#'         intensity = c(100, 80.1)
#'     ),
#'     data.frame(
#'         rtime = c(12.4, 12.8, 13.2, 14.6),
#'         intensity = c(123.3, 153.6, 2354.3, 243.4)
#'     )
#' )
#'
#' cbm <- backendInitialize(cbm, chromData = cdata, peaksData = pdata)
#' cbm
#'
NULL

#' @noRd
setClass("ChromBackendMemory",
    contains = "ChromBackend",
    slots = c(
        chromData = "data.frame",
        peaksData = "list"
    ),
    prototype = prototype(
        chromData = fillCoreChromVariables(data.frame()),
        peaksData = list(.EMPTY_PEAKS_DATA),
        readonly = FALSE,
        version = "0.1"
    )
)

#' @rdname ChromBackendMemory
#'
#' @importFrom methods new
#' @export ChromBackendMemory
ChromBackendMemory <- function() {
    new("ChromBackendMemory")
}

#' @rdname ChromBackendMemory
setMethod(
    "backendInitialize", "ChromBackendMemory",
    function(object,
    chromData = fillCoreChromVariables(data.frame()),
    peaksData = list(.EMPTY_PEAKS_DATA), ...) {
        if (!is(chromData, "data.frame")) {
            stop(
                "'chromData' needs to be a 'data.frame' with the general",
                " chromatogram variables"
            )
        }
        n_cd <- nrow(chromData)
        if (n_cd) {
            if (is.null(chromData$dataOrigin)) {
                chromData$dataOrigin <- NA_character_
            }
            validChromData(chromData)
            chromData <- chromData[, !vapply(chromData, function(x) {
                all(is.na(x))
            }, logical(1)), drop = FALSE]
        } else {
            chromData <- fillCoreChromVariables(data.frame())
        }
        object@chromData <- chromData
        if (length(peaksData) > 1 || length(peaksData) == n_cd) {
            validPeaksData(peaksData)
        } else {
            peaksData <- replicate(n_cd, .EMPTY_PEAKS_DATA,
                simplify = FALSE
            )
        }
        object@peaksData <- peaksData
        object
    }
)

#' @rdname hidden_aliases
setMethod("backendMerge", "ChromBackendMemory", function(object, ...) {
    object <- unname(c(object, ...))
    not_empty <- lengths(object) > 0
    if (any(not_empty)) {
        res <- .df_combine(object[not_empty])
    } else {
        res <- object[[1L]]
    }
    validChromData(res@chromData)
    validPeaksData(res@peaksData)
    res
})

#' @rdname hidden_aliases
#' @description This method returns the chromatographic data stored in the
#' backend. If not specified otherwise it will return all defined columns in
#' the chromData slot as well as adding the `coreChromVariables` missing with
#' NA values.
setMethod(
    "chromData", "ChromBackendMemory",
    function(object, columns = chromVariables(object), drop = FALSE) {
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

#' @rdname hidden_aliases
setReplaceMethod("chromData", "ChromBackendMemory", function(object, value) {
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

#' @rdname hidden_aliases
setMethod("chromVariables", "ChromBackendMemory", function(object) {
    union(names(object@chromData), names(coreChromVariables()))
})

#' @rdname hidden_aliases
setMethod(
    "peaksData", "ChromBackendMemory",
    function(object, columns = peaksVariables(object),
    drop = FALSE, ...) {
        if (!any(peaksVariables(object) %in% columns)) {
            stop(
                "Some of the requested peaks variables are not",
                " available"
            )
        }
        if (identical(
            as.vector(peaksVariables(object)),
            as.vector(columns)
        )) {
            return(object@peaksData)
        }
        lapply(object@peaksData, function(x) x[, columns, drop = drop])
    }
)

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "ChromBackendMemory", function(object, value) {
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

#' @rdname hidden_aliases
setMethod("peaksVariables", "ChromBackendMemory", function(object) {
    union(names(object@peaksData[[1]]), names(corePeaksVariables()))
})

#' @rdname hidden_aliases
#' @export
setMethod("isReadOnly", "ChromBackendMemory", function(object) FALSE)

#' @importFrom utils capture.output head
#' @rdname hidden_aliases
setMethod("show", "ChromBackendMemory", function(object) {
    cpd <- object@chromData
    cat(class(object), "with", nrow(cpd), "chromatograms\n")
    if (nrow(cpd)) {
        cpd <- fillCoreChromVariables(cpd)
        cpd <- cpd[, c("chromIndex", "msLevel", "mz"), drop = FALSE]
        txt <- capture.output(print(head(cpd)))
        cat(txt, sep = "\n")
        cp_cols <- object@chromData[
            ,
            !(colnames(object@chromData) %in%
                colnames(cpd))
        ]
        cat("...", length(cp_cols), "more  chromatogram variables/columns\n")
        cat("...", ncol(object@peaksData[[1]]), "peaksData variables\n")
    }
})

#' @rdname hidden_aliases
#' @export
setMethod(
    "supportsSetBackend", "ChromBackendMemory",
    function(object, ...) TRUE
)

#' @importFrom MsCoreUtils i2index
#' @importMethodsFrom S4Vectors [ [<-
#' @rdname hidden_aliases
setMethod("[", "ChromBackendMemory", function(x, i, j, ..., drop = FALSE) {
    if (!length(i)) {
        return(ChromBackendMemory())
    }
    i <- i2index(i, length = length(x))
    x@chromData <- x@chromData[i, ]
    x@peaksData <- x@peaksData[i]
    x
})

#' @rdname hidden_aliases
setMethod("$", "ChromBackendMemory", function(x, name) {
    if (name %in% union(chromVariables(x), names(coreChromVariables()))) {
        res <- chromData(x, columns = name, drop = TRUE)
    } else if (name %in% peaksVariables(x)) {
        res <- peaksData(x, columns = name, drop = TRUE)
    } else {
        stop("The requested variable '", name, "' is not available")
    }
    res
})

#' @rdname hidden_aliases
setReplaceMethod("$", "ChromBackendMemory", function(x, name, value) {
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
        for (i in seq_along(value)) peaksData(x)[[i]][[name]] <- value[[i]]
    } else {
        chromData(x)[, name] <- value
    }
    x
})
