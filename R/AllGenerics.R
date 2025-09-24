setGeneric("Chromatograms", function(object, ...) {
    standardGeneric("Chromatograms")
})
setGeneric("chromData", function(object, ...) standardGeneric("chromData"))
setGeneric("chromData<-", function(object, value) {
    standardGeneric("chromData<-")
})
setGeneric("chromExtract", function(object,peak_table, by, ...) standardGeneric("chromExtract"))
setGeneric("chromIndex", function(object, ...) standardGeneric("chromIndex"))
setGeneric("chromIndex<-", function(object, value) {
    standardGeneric("chromIndex<-")
})
setGeneric("chromVariables", function(object, ...) {
    standardGeneric("chromVariables")
})
setGeneric("factorize", function(object, ...) standardGeneric("factorize"))
setGeneric("filterChromData", function(object, ...) {
    standardGeneric("filterChromData")
})
setGeneric("filterPeaksData", function(object, ...) {
    standardGeneric("filterPeaksData")
})
setGeneric("mzMax", function(object, ...) standardGeneric("mzMax"))
setGeneric("mzMax<-", function(object, value) standardGeneric("mzMax<-"))
setGeneric("mzMin", function(object, ...) standardGeneric("mzMin"))
setGeneric("mzMin<-", function(object, value) standardGeneric("mzMin<-"))
setGeneric("precursorMzMin", function(object, ...) {
    standardGeneric("precursorMzMin")
})
setGeneric("precursorMzMin<-", function(object, value) {
    standardGeneric("precursorMzMin<-")
})
setGeneric("precursorMzMax", function(object, ...) {
    standardGeneric("precursorMzMax")
})
setGeneric("precursorMzMax<-", function(object, value) {
    standardGeneric("precursorMzMax<-")
})
setGeneric("productMzMax", function(object, ...) {
    standardGeneric("productMzMax")
})
setGeneric("productMzMax<-", function(object, value) {
    standardGeneric("productMzMax<-")
})
setGeneric("productMzMin", function(object, ...) {
    standardGeneric("productMzMin")
})
setGeneric("productMzMin<-", function(object, value) {
    standardGeneric("productMzMin<-")
})
setGeneric("imputePeaksData", function(object, ...)
    standardGeneric("imputePeaksData"))
#' @rdname hidden_aliases
setGeneric("reset", function(object, ...) {
    standardGeneric("reset")
}) ## needs to be moved to ProtGenerics
