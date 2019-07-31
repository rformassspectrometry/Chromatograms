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
setGeneric("isReadOnly", function(object, ...)
    standardGeneric("isReadOnly"))
