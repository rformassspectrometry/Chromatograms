#' @include hidden_aliases.R
NULL

#' @rdname hidden_aliases
setGeneric("backendInitialize",
           def = function(object, ...)
               standardGeneric("backendInitialize"),
           valueClass = "ChromBackend"
           )
#' @rdname hidden_aliases
setGeneric("backendMerge", def = function(object, ...)
    standardGeneric("backendMerge"),
    valueClass = "ChromBackend")
#' @rdname hidden_aliases
setGeneric("chromData", function(object, ...) standardGeneric("chromData"))
#' @rdname hidden_aliases
setGeneric("chromData<-", function(object, value)
    standardGeneric("chromData<-"))
#' @rdname hidden_aliases
setGeneric("chromNames", function(object, ...) standardGeneric("chromNames"))
#' @rdname hidden_aliases
setGeneric("chromNames<-", function(object, value)
    standardGeneric("chromNames<-"))
#' @rdname hidden_aliases
setGeneric("chromVariables", function(object, ...)
    standardGeneric("chromVariables"))
#' @rdname hidden_aliases
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
#' @rdname hidden_aliases
setGeneric("mzMax", function(object, ...) standardGeneric("mzMax"))
#' @rdname hidden_aliases
setGeneric("mzMax<-", function(object, value) standardGeneric("mzMax<-"))
#' @rdname hidden_aliases
setGeneric("mzMin", function(object, ...) standardGeneric("mzMin"))
#' @rdname hidden_aliases
setGeneric("mzMin<-", function(object, value) standardGeneric("mzMin<-"))
#' @rdname hidden_aliases
setGeneric("precursorMzMin", function(object, ...)
    standardGeneric("precursorMzMin"))
#' @rdname hidden_aliases
setGeneric("precursorMzMin<-", function(object, value)
    standardGeneric("precursorMzMin<-"))
#' @rdname hidden_aliases
setGeneric("precursorMzMax", function(object, ...)
    standardGeneric("precursorMzMax"))
#' @rdname hidden_aliases
setGeneric("precursorMzMax<-", function(object, value)
    standardGeneric("precursorMzMax<-"))
#' @rdname hidden_aliases
setGeneric("productMzMax", function(object, ...)
    standardGeneric("productMzMax"))
#' @rdname hidden_aliases
setGeneric("productMzMax<-", function(object, value)
    standardGeneric("productMzMax<-"))
#' @rdname hidden_aliases
setGeneric("productMzMin", function(object, ...)
    standardGeneric("productMzMin"))
#' @rdname hidden_aliases
setGeneric("productMzMin<-", function(object, value)
    standardGeneric("productMzMin<-"))
#' @rdname hidden_aliases
setGeneric("selectChromVariables", function(object, ...)
    standardGeneric("selectChromVariables"))
