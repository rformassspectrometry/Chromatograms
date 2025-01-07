setGeneric("chromData", function(object, ...) standardGeneric("chromData"))
#' @rdname hidden_aliases
setGeneric("chromData<-", function(object, value)
    standardGeneric("chromData<-"))
setGeneric("chromIndex", function(object, ...) standardGeneric("chromIndex"))
#' @rdname hidden_aliases
setGeneric("chromIndex<-", function(object, value)
    standardGeneric("chromIndex<-"))
setGeneric("chromVariables", function(object, ...)
    standardGeneric("chromVariables"))
setGeneric("filterChromData", function(object, ...)
    standardGeneric("filterChromData"))
setGeneric("mzMax", function(object, ...) standardGeneric("mzMax"))
#' @rdname hidden_aliases
setGeneric("mzMax<-", function(object, value) standardGeneric("mzMax<-"))
setGeneric("mzMin", function(object, ...) standardGeneric("mzMin"))
#' @rdname hidden_aliases
setGeneric("mzMin<-", function(object, value) standardGeneric("mzMin<-"))
setGeneric("precursorMzMin", function(object, ...)
    standardGeneric("precursorMzMin"))
#' @rdname hidden_aliases
setGeneric("precursorMzMin<-", function(object, value)
    standardGeneric("precursorMzMin<-"))
setGeneric("precursorMzMax", function(object, ...)
    standardGeneric("precursorMzMax"))
#' @rdname hidden_aliases
setGeneric("precursorMzMax<-", function(object, value)
    standardGeneric("precursorMzMax<-"))
setGeneric("productMzMax", function(object, ...)
    standardGeneric("productMzMax"))
#' @rdname hidden_aliases
setGeneric("productMzMax<-", function(object, value)
    standardGeneric("productMzMax<-"))
setGeneric("productMzMin", function(object, ...)
    standardGeneric("productMzMin"))
#' @rdname hidden_aliases
setGeneric("productMzMin<-", function(object, value)
    standardGeneric("productMzMin<-"))
#' @rdname hidden_aliases
setGeneric("reset", function(object, ...)
    standardGeneric("reset")) ## needs to be moved to ProtGenerics
