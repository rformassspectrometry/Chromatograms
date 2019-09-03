#' @include hidden_aliases.R
NULL

#' @noRd
#'
#' @importClassesFrom S4Vectors DataFrame Rle
#'
#' @importFrom S4Vectors DataFrame Rle
#'
#' @author Johannes Rainer
#'
#' @exportClass ChromBackendDataFrame
setClass("ChromBackendDataFrame",
         contains = "ChromBackend", slots = c(chromData = "DataFrame"),
         prototype = prototype(chromData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

setValidity("ChromBackendDataFrame", function(object) {
    msg <- .valid_chrom_data_required_columns(object@chromData)
    if (length(msg))
        return(msg)
    msg <- c(
        .valid_column_datatype(object@chromData, .CHROMATOGRAMS_DATA_COLUMNS),
        .valid_intensity_column(object@chromData),
        .valid_rtime_column(object@chromData),
        .valid_intensity_rtime_columns(object@chromData))
    if (is.null(msg)) TRUE
    else msg
})

#' @rdname hidden_aliases
#'
#' @importFrom methods show
#'
#' @importFrom utils capture.output
setMethod("show", "ChromBackendDataFrame", function(object) {
    chd <- chromData(object, c("msLevel", "mz", "chromIndex"))
    cat(class(object), "with", nrow(chd), "chromatograms\n")
    if (nrow(chd)) {
        txt <- capture.output(print(chd))
        cat(txt[-1], sep = "\n")
        ch_cols <- chromVariables(object)
        cat(" ...", length(ch_cols) - 3, "more variables/columns.\n")
    }
})

#' @importMethodsFrom S4Vectors $ $<-
#'
#' @importClassesFrom IRanges NumericList
#'
#' @importFrom MsCoreUtils asRleDataFrame
#'
#' @importFrom IRanges NumericList
#'
#' @rdname hidden_aliases
setMethod("backendInitialize", signature = "ChromBackendDataFrame",
          function(object, chromData, ...) {
              if (missing(chromData)) chromData <- DataFrame()
              if (is.data.frame(chromData))
                  chromData <- DataFrame(chromData)
              if (!is(chromData, "DataFrame"))
                  stop("'chromData' has to be a 'DataFrame'")
              if (!nrow(chromData))
                  return(object)
              chromData$dataStorage <- "<memory>"
              chromData <- asRleDataFrame(
                  chromData, columns = c("dataStorage", "dataOrigin"))
              if (nrow(chromData) && !is(chromData$rtime, "NumericList"))
                  chromData$rtime <- NumericList(chromData$rtime,
                                                 compress = FALSE)
              if (nrow(chromData) && !is(chromData$intensity, "NumericList"))
                  chromData$intensity <- NumericList(chromData$intensity,
                                                     compress = FALSE)
              object@chromData <- chromData
              validObject(object)
              object
          })

#' @rdname hidden_aliases
setMethod("backendMerge", "ChromBackendDataFrame", function(object, ...) {
    object <- unname(c(object, ...))
    object <- object[lengths(object) > 0]
    res <- .combine_chrom_backend_data_frame(object)
    validObject(res)
    res
})

#' @rdname hidden_aliases
#'
#' @importFrom methods as
#'
#' @importFrom S4Vectors SimpleList
#'
#' @importFrom MsCoreUtils asVectorDataFrame
#'
#' @importMethodsFrom S4Vectors lapply
setMethod("chromData", "ChromBackendDataFrame",
          function(object, columns = chromVariables(object)) {
              df_columns <- intersect(columns, colnames(object@chromData))
              res <- object@chromData[, df_columns, drop = FALSE]
              other_columns <- setdiff(columns, colnames(object@chromData))
              if (length(other_columns)) {
                  other_res <- lapply(other_columns, .get_rle_column,
                                      x = object@chromData)
                  names(other_res) <- other_columns
                  is_rtime_int <- names(other_res) %in% c("rtime", "intensity")
                  if (!all(is_rtime_int))
                      res <- cbind(res, as(other_res[!is_rtime_int], "DataFrame"))
                  if (any(names(other_res) == "rtime"))
                      res$rtime <- if (length(other_res$rtime)) other_res$rtime
                                else NumericList(compress = FALSE)
                  if (any(names(other_res) == "intensity"))
                      res$intensity <- if (length(other_res$intensity)) other_res$intensity
                                       else NumericList(compress = FALSE)
              }
              asVectorDataFrame(res[, columns, drop = FALSE])
          })

#' @rdname hidden_aliases
#'
#' @importFrom MsCoreUtils asRleDataFrame
setReplaceMethod("chromData", "ChromBackendDataFrame", function(object, value) {
    if (inherits(value, "DataFrame")) {
        if (length(object) && nrow(value) != length(object))
            stop("'value' has to be a 'DataFrame' with ", length(object), " rows.")
        if (!is(value$rtime, "NumericList"))
            value$rtime <- NumericList(value$rtime, compress = FALSE)
        if (!is(value$intensity, "NumericList"))
            value$intensity <- NumericList(value$intensity, compress = FALSE)
        if (is.null(value$dataStorage))
            value$dataStorage <- "<memory>"
    } else {
        if (length(value) == 1)
            value <- rep.int(value, length(object))
        if (length(value) != length(object))
            stop("length of 'value' has to be ", length(object))
    }
    object@chromData <- asRleDataFrame(value, columns = c("dataStorage",
                                                          "dataOrigin"))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("chromIndex", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "chromIndex")
})

#' @rdname hidden_aliases
setMethod("chromNames", "ChromBackendDataFrame", function(object) {
    rownames(object@chromData)
})

#' @rdname hidden_aliases
setReplaceMethod("chromNames", "ChromBackendDataFrame", function(object, value) {
    rownames(object@chromData) <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("chromVariables", "ChromBackendDataFrame", function(object) {
    unique(c(names(.CHROMATOGRAMS_DATA_COLUMNS), colnames(object@chromData)))
})

#' @rdname hidden_aliases
setMethod("collisionEnergy", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "collisionEnergy")
})

#' @rdname hidden_aliases
#'
#' @importFrom MsCoreUtils asRle
setReplaceMethod("collisionEnergy", "ChromBackendDataFrame", function(object,
                                                                      value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$collisionEnergy <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataOrigin", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "dataOrigin")
})

#' @rdname hidden_aliases
#'
setReplaceMethod("dataOrigin", "ChromBackendDataFrame", function(object,
                                                                 value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@chromData$dataOrigin <- asRle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("dataStorage", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "dataStorage")
})

#' @rdname hidden_aliases
setReplaceMethod("dataStorage", "ChromBackendDataFrame", function(object,
                                                                  value) {
    if (!is.character(value) || length(value) != length(object))
        stop("'value' has to be a 'character' of length ", length(object))
    object@chromData$dataStorage <- asRle(as.character(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("intensity", "ChromBackendDataFrame", function(object) {
    if (any(colnames(object@chromData) == "intensity"))
        object@chromData$intensity
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep.int(1, length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("intensity", "ChromBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == lengths(object$rtime)))
        stop("lengths of 'value' has to match the number of data pairs ",
             "(i.e. pairsCount(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@chromData$intensity <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("isEmpty", "ChromBackendDataFrame", function(x) {
    lengths(intensity(x)) == 0
})

#' @rdname hidden_aliases
setMethod("length", "ChromBackendDataFrame", function(x) {
    nrow(x@chromData)
})

#' @rdname hidden_aliases
setMethod("msLevel", "ChromBackendDataFrame", function(object, ...) {
    .get_rle_column(object@chromData, "msLevel")
})

#' @rdname hidden_aliases
#'
#' @importFrom MsCoreUtils asInteger
setReplaceMethod("msLevel", "ChromBackendDataFrame", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$msLevel <- asRle(asInteger(value))
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("mz", "ChromBackendDataFrame", function(object, ...) {
    .get_rle_column(object@chromData, "mz")
})

#' @rdname hidden_aliases
setReplaceMethod("mz", "ChromBackendDataFrame", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$mz <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("mzMax", "ChromBackendDataFrame", function(object, ...) {
    vals <- .get_rle_column(object@chromData, "mzMax")
    if (all(is.na(vals)))
        mz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("mzMax", "ChromBackendDataFrame", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$mzMax <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("mzMin", "ChromBackendDataFrame", function(object, ...) {
    vals <- .get_rle_column(object@chromData, "mzMin")
    if (all(is.na(vals)))
        mz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("mzMin", "ChromBackendDataFrame", function(object, value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$mzMin <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("pairs", "ChromBackendDataFrame", function(x, ...) {
    mapply(cbind, rtime = rtime(x), intensity = intensity(x),
           SIMPLIFY = FALSE, USE.NAMES = FALSE)
})

#' @rdname hidden_aliases
setReplaceMethod("pairs", "ChromBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "SimpleList")))
        stop("'value' has to be a list-like object")
    if (length(value) != length(object))
        stop("Length of 'value' has to match length of 'object'")
    vals <- lapply(value, "[", , 1L)
    if (!is(vals, "NumericList"))
        vals <- NumericList(vals, compress = FALSE)
    object@chromData$rtime <- vals
    vals <- lapply(value, "[", , 2L)
    if (!is(vals, "NumericList"))
        vals <- NumericList(vals, compress = FALSE)
    object@chromData$intensity <- vals
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("pairsCount", "ChromBackendDataFrame", function(object) {
    lengths(intensity(object))
})

#' @rdname hidden_aliases
setMethod("precursorMz", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "precursorMz")
})

#' @rdname hidden_aliases
setReplaceMethod("precursorMz", "ChromBackendDataFrame", function(object,
                                                                  value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$precursorMz <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precursorMzMax", "ChromBackendDataFrame", function(object) {
    vals <- .get_rle_column(object@chromData, "precursorMzMax")
    if (all(is.na(vals)))
        precursorMz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("precursorMzMax", "ChromBackendDataFrame", function(object,
                                                                  value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$precursorMzMax <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("precursorMzMin", "ChromBackendDataFrame", function(object) {
    vals <- .get_rle_column(object@chromData, "precursorMzMin")
    if (all(is.na(vals)))
        precursorMz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("precursorMzMin", "ChromBackendDataFrame", function(object,
                                                                     value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$precursorMzMin <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("productMz", "ChromBackendDataFrame", function(object) {
    .get_rle_column(object@chromData, "productMz")
})

#' @rdname hidden_aliases
setReplaceMethod("productMz", "ChromBackendDataFrame", function(object,
                                                                value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$productMz <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("productMzMax", "ChromBackendDataFrame", function(object) {
    vals <- .get_rle_column(object@chromData, "productMzMax")
    if (all(is.na(vals)))
        productMz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("productMzMax", "ChromBackendDataFrame", function(object,
                                                                   value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$productMzMax <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("productMzMin", "ChromBackendDataFrame", function(object) {
    vals <- .get_rle_column(object@chromData, "productMzMin")
    if (all(is.na(vals)))
        productMz(object)
    else vals
})

#' @rdname hidden_aliases
setReplaceMethod("productMzMin", "ChromBackendDataFrame", function(object,
                                                                  value) {
    if (!is.numeric(value) || length(value) != length(object))
        stop("'value' has to be a 'numeric' of length ", length(object))
    object@chromData$productMzMin <- asRle(value)
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("rtime", "ChromBackendDataFrame", function(object) {
    if (any(colnames(object@chromData) == "rtime"))
        object@chromData$rtime
    else {
        lst <- NumericList(numeric(), compress = FALSE)
        lst[rep.int(1, length(object))]
    }
})

#' @rdname hidden_aliases
setReplaceMethod("rtime", "ChromBackendDataFrame", function(object, value) {
    if (!(is.list(value) || inherits(value, "NumericList")))
        stop("'value' has to be a list or NumericList")
    if (length(value) != length(object))
        stop("length of 'value' has to match the length of 'object'")
    if (!all(lengths(value) == pairsCount(object)))
        stop("lengths of 'value' has to match the number of data points ",
             "(i.e. pairsCount(object))")
    if (!is(value, "NumericList"))
        value <- NumericList(value, compress = FALSE)
    object@chromData$rtime <- value
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod("selectChromVariables", "ChromBackendDataFrame",
          function(object, chromVariables = chromVariables(object)) {
              if (!all(chromVariables %in% chromVariables(object)))
                  stop("Chromarogram variables ",
                       paste(chromVariables[!(chromVariables %in%
                                              chromVariables(object))],
                             collapse = ", "), " not available")
              keep <- chromVariables[chromVariables %in%
                                     colnames(object@chromData)]
              if (length(keep))
                  object@chromData <- object@chromData[, keep, drop = FALSE]
              msg <- .valid_chrom_data_required_columns(object@chromData)
              if (length(msg))
                  stop(msg)
              validObject(object)
              object
})

#' @rdname hidden_aliases
setMethod("$", "ChromBackendDataFrame", function(x, name) {
    .get_rle_column(x@chromData, column = name)
})

#' @rdname hidden_aliases
setReplaceMethod("$", "ChromBackendDataFrame", function(x, name, value) {
    if (is.list(value) && any(c("rtime", "intensity") == name))
        value <- NumericList(value, compress = FALSE)
    if (name == "dataStorage")
        value <- Rle(value)
    x@chromData[[name]] <- asRle(value)
    validObject(x)
    x
})

#### ---------------------------------------------------------------------------
##
##                      FILTERING AND SUBSETTING
##
#### ---------------------------------------------------------------------------

#' @importMethodsFrom S4Vectors [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @rdname hidden_aliases
setMethod("[", "ChromBackendDataFrame", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j))
        stop("Subsetting by column ('j = ", j, "' is not supported")
    i <- i2index(i, length(x), rownames(x@chromData))
    x@chromData <- x@chromData[i, , drop = FALSE]
    validObject(x)
    x
})

#' @rdname hidden_aliases
setMethod("filterDataOrigin", "ChromBackendDataFrame",
          function(object, dataOrigin = character(), ...) {
              if (length(dataOrigin)) {
                  object <- object[dataOrigin(object) %in% dataOrigin]
                  if (is.unsorted(dataOrigin))
                      object[order(match(dataOrigin(object), dataOrigin))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterDataStorage", "ChromBackendDataFrame",
          function(object, dataStorage = character(), ...) {
              if (length(dataStorage)) {
                  object <- object[dataStorage(object) %in% dataStorage]
                  if (is.unsorted(dataStorage))
                      object[order(match(dataStorage(object), dataStorage))]
                  else object
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterMsLevel", "ChromBackendDataFrame",
          function(object, msLevel = integer(), ...) {
              if (length(msLevel))
                  object[msLevel(object) %in% msLevel]
              else object
          })

#' @rdname hidden_aliases
setMethod("filterMz", "ChromBackendDataFrame",
          function(object, mz = numeric(), ...) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(mzMax(object) >= mz[1] &
                                mzMin(object) <= mz[2])
                  object[keep]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterPrecursorMz", "ChromBackendDataFrame",
          function(object, mz = numeric(), ...) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(precursorMzMax(object) >= mz[1] &
                                precursorMzMin(object) <= mz[2])
                  object[keep]
              } else object
          })

#' @rdname hidden_aliases
setMethod("filterProductMz", "ChromBackendDataFrame",
          function(object, mz = numeric(), ...) {
              if (length(mz)) {
                  mz <- range(mz)
                  keep <- which(productMzMax(object) >= mz[1] &
                                productMzMin(object) <= mz[2])
                  object[keep]
              } else object
})
