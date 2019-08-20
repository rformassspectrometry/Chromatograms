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
setGeneric("dataOrigin", function(object, ...) standardGeneric("dataOrigin"))
#' @rdname hidden_aliases
setGeneric("dataOrigin<-", function(object, value)
    standardGeneric("dataOrigin<-"))
#' @rdname hidden_aliases
setGeneric("dataStorage", function(object, ...) standardGeneric("dataStorage"))
#' @rdname hidden_aliases
setGeneric("dataStorage<-", function(object, value)
    standardGeneric("dataStorage<-"))
#' @rdname hidden_aliases
setGeneric("filterDataOrigin", function(object, ...)
    standardGeneric("filterDataOrigin"))
#' @rdname hidden_aliases
setGeneric("filterDataStorage", function(object, ...)
    standardGeneric("filterDataStorage"))
#' @rdname hidden_aliases
setGeneric("filterMsLevel", function(object, ...)
    standardGeneric("filterMsLevel"))
#' @rdname hidden_aliases
setGeneric("filterMz", function(object, ...) standardGeneric("filterMz"))
#' @rdname hidden_aliases
setGeneric("filterPolarity", function(object, ...)
    standardGeneric("filterPolarity"))
#' @rdname hidden_aliases
setGeneric("filterPrecursorMz", function(object, ...)
    standardGeneric("filterPrecursorMz"))
#' @rdname hidden_aliases
setGeneric("filterProductMz", function(object, ...)
    standardGeneric("filterProductMz"))
#' @rdname hidden_aliases
setGeneric("filterRt", function(object, ...) standardGeneric("filterRt"))
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
setGeneric("pairs", function(object, ...) standardGeneric("pairs"))
#' @rdname hidden_aliases
setGeneric("pairs<-", function(object, value) standardGeneric("pairs<-"))
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
setGeneric("productMz", function(object, ...) standardGeneric("productMz"))
#' @rdname hidden_aliases
setGeneric("productMz<-", function(object, value) standardGeneric("productMz<-"))
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
