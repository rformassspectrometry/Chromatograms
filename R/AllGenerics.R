#' @rdname ChromBackend
setGeneric("chromData", function(object, ...) standardGeneric("chromData"))
#' @rdname ChromBackend
setGeneric("chromData<-", function(object, value)
    standardGeneric("chromData<-"))
#' @rdname ChromBackend
setGeneric("chromIndex", function(object, ...) standardGeneric("chromIndex"))
#' @rdname ChromBackend
setGeneric("chromIndex<-", function(object, value)
    standardGeneric("chromIndex<-"))
#' @rdname ChromBackend
setGeneric("chromVariables", function(object, ...)
    standardGeneric("chromVariables"))
#' @rdname ChromBackend
setGeneric("mzMax", function(object, ...) standardGeneric("mzMax"))
#' @rdname ChromBackend
setGeneric("mzMax<-", function(object, value) standardGeneric("mzMax<-"))
#' @rdname ChromBackend
setGeneric("mzMin", function(object, ...) standardGeneric("mzMin"))
#' @rdname ChromBackend
setGeneric("mzMin<-", function(object, value) standardGeneric("mzMin<-"))
#' @rdname ChromBackend
setGeneric("precursorMzMin", function(object, ...)
    standardGeneric("precursorMzMin"))
#' @rdname ChromBackend
setGeneric("precursorMzMin<-", function(object, value)
    standardGeneric("precursorMzMin<-"))
#' @rdname ChromBackend
setGeneric("precursorMzMax", function(object, ...)
    standardGeneric("precursorMzMax"))
#' @rdname ChromBackend
setGeneric("precursorMzMax<-", function(object, value)
    standardGeneric("precursorMzMax<-"))
#' @rdname ChromBackend
setGeneric("productMzMax", function(object, ...)
    standardGeneric("productMzMax"))
#' @rdname ChromBackend
setGeneric("productMzMax<-", function(object, value)
    standardGeneric("productMzMax<-"))
#' @rdname ChromBackend
setGeneric("productMzMin", function(object, ...)
    standardGeneric("productMzMin"))
#' @rdname ChromBackend
setGeneric("productMzMin<-", function(object, value)
    standardGeneric("productMzMin<-"))
#' @rdname ChromBackend
setGeneric("reset", function(object, ...)
    standardGeneric("reset")) ## needs to be moved to ProtGenerics
