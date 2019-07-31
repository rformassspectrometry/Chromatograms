#' @include hidden_aliases.R
NULL

#' @title Chromatographic data backends
#'
#' @description
#'
#' `ChromBackend` is a virtual class that defines what different backends need
#' to provide to be used by the `Chromatograms` package and classes.
#'
#' @name ChromBackend
#'
#' @author Johannes Rainer, Sebastian Gibb, Laurent Gatto
#'
#' @md
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

## chromData, chromData<-
## chromNames, chromNames<-
## chromVariables

#' @exportMethod dataOrigin
#'
#' @rdname ChromBackend
setMethod("dataOrigin", "ChromBackend", function(object) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataOrigin<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataOrigin", "ChromBackend", function(object, value) {
    stop("Not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage
#'
#' @rdname ChromBackend
setMethod("dataStorage", "ChromBackend", function(object) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
})

#' @exportMethod dataStorage<-
#'
#' @rdname ChromBackend
setReplaceMethod("dataStorage", "ChromBackend", function(object, value) {
    stop("Method 'dataStorage' is not implemented for ", class(object), ".")
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

## mzmin, mzmin<- ... precursorMz? productMz?
## mzmax, mzmax<-
## pairs, pairs<-
## rtime, rtime<-
## [
## $, $<-
