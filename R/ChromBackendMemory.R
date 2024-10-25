#' @include helpers.R
#' @include hidden_aliases.R
#' @include ChromBackend.R
NULL

#' @title Improved in-memory Chromatographic data backend
#'
#' @description
#' This backend is used to store chromatographic data in memory. It is
#' particularly useful for small datasets or for testing purposes.
#'
#'
#' @param chromData For `backendInitialize()`, a `data.frame` with the
#'     chromatographic data. If not provided (or if empty), a default
#'     `data.frame` with the core chromatographic variables will be created.
#'
#' @param peaksData For `backendInitialize()`, a `list` of `data.frame` with
#'     the peaks data. If not provided (or if empty), a default `list` of empty
#'     `data.frame` with the core peaks variables will be created. The length
#'     of the list should match the number of chromatograms in the `chromData`
#'     parameter.
#'
#'
#' @author Philippine Louail
#'
#' @exportClass ChromBackendMemory
#'
#' @noRd
NULL

#' ## Some definitions
#' - chromData: `data.frame` with arbitrary spectra annotations.
#' - peaksData: `list` of `data.frame` with retention time and intensity
#' values.
#'
#' @noRd
setClass("ChromBackendMemory",
         contains = "ChromBackend",
         slots = c(chromData = "data.frame",
                   peaksData = "list"),
         prototype = prototype(chromData = fillCoreChromVariables(data.frame()),
                               peaksData = list(.EMPTY_PEAKS_DATA),
                               readonly = FALSE,
                               version = "0.1"))

#' @rdname ChromBackend
setMethod("backendInitialize", "ChromBackendMemory",
          function(object,
                   chromData = fillCoreChromVariables(data.frame()),
                   peaksData = list(.EMPTY_PEAKS_DATA)) {

            # Ensure chromData is a data.frame
            if (!is(chromData, "data.frame"))
              stop("'chromData' needs to be a 'data.frame' with the general",
                   " chromatogram variables")

            n_cd <- nrow(chromData)
            # If chromData is not empty, fill default values if needed
            if (n_cd) {
              if (is.null(chromData$dataStorage))
                chromData$dataStorage <- "<memory>"
              if (is.null(chromData$dataOrigin))
                chromData$dataOrigin <- "<user provided>"

              # Validate and clean up chromData
              validChromData(chromData)
              chromData <- chromData[, !vapply(chromData,
                                                     function(x) all(is.na(x)),
                                                     logical(1)), drop = FALSE]
            } else {
              chromData <- fillCoreChromVariables(data.frame())
            }

            # Assign chromData to the object@chromData slot
            object@chromData <- chromData

            # Handle peaksData
            if (length(peaksData) > 1) {
              validPeaksData(peaksData)
            } else {
              peaksData <- replicate(n_cd, .EMPTY_PEAKS_DATA, simplify = FALSE)
            }
            # Assign peaksData to the object@peaksData slot
            object@peaksData <- peaksData

            object
          })

#' @rdname hidden_aliases
setMethod("backendMerge", "ChromBackendMemory", function(object, ...) {
  object <- unname(c(object, ...))
  not_empty <- lengths(object) > 0
  if (any(not_empty))
    res <- .df_combine(object[not_empty])
  else res <- object[[1L]]
  validChromData(res@chromData)
  validPeaksData(res@peaksData)
  res
})

#' @rdname hidden_aliases
#' @description This method returns the chromatographic data stored in the
#' backend. If not specified otherwise it will return all defined column in the
#' chromData slot as well as dding the coreChromVariables missing with NA values.
#'
#' @importMethodsFrom S4Vectors sapply
setMethod("chromData", "ChromBackendMemory",
          function(object,
                   columns = union(names(coreChromVariables()),
                                   chromVariables(object)),
                   drop = FALSE) {
            if (!any(chromVariables(object) %in% columns) &&
                !any(names(coreChromVariables()) %in% columns))
              stop("Some of the requested Chromatogram variables are not ",
                   "available")
            if (all(columns %in% chromVariables(object)))
              res <- object@chromData[, columns, drop = drop]
            else {
              res <- fillCoreChromVariables(object@chromData)
              res <- res[, columns, drop = drop]
            }
            res
          })

#' @rdname hidden_aliases
setReplaceMethod("chromData", "ChromBackendMemory", function(object, value) {
  if (!inherits(value, "data.frame"))
    stop("'value' is expected to be a 'data.frame'")
  if (length(object) && length(object) != nrow(value))
    stop("'value' has to be a 'data.frame' with ", length(object), " rows")
  validChromData(value)
  object@chromData <- value
  object
})

#' @rdname hidden_aliases
setMethod("chromVariables", "ChromBackendMemory", function(object) {
  names(object@chromData)
})

#' @rdname hidden_aliases
setMethod("isReadOnly", "ChromBackendMemory", function(object) FALSE)

#' @rdname hidden_aliases
setMethod("peaksData", "ChromBackendMemory",
          function(object, columns = peaksVariables(object), drop = FALSE) {
            if (!all(columns %in% peaksVariables(object)))
              stop("Some of the requested peaks variables are not available")
            res <- lapply(object@peaksData, function(x) x[, columns,
                                                          drop = drop])
            res
          })

#' @rdname hidden_aliases
setReplaceMethod("peaksData", "ChromBackendMemory", function(object, value) {
  if (!is.list(value))
    stop("'value' is expected to be a list")
  if (length(object) && length(object) != length(value))
    stop("'value' has to be a list with ", length(object), " elements")
  validPeaksData(value)
  object@peaksData <- value
  object
})

#' @rdname hidden_aliases
setMethod("peaksVariables", "ChromBackendMemory", function(object) {
  names(object@peaksData[[1]])
})

#' @importFrom utils capture.output
#' @rdname hidden_aliases
setMethod("show", "ChromBackendMemory", function(object){
  cpd <- object@chromData
  cat(class(object), "with", nrow(cpd), "chromatograms\n")
  if (nrow(cpd)) {
    cpd <- cpd[, c("chromIndex", "msLevel", "mz")]
    txt <- capture.output(print(cpd))
    cat(txt, sep = "\n")
    cp_cols <- object@chromData[,
                                !(colnames(object@chromData) %in% colnames(cpd))]
    cat("...", length(cp_cols), "more variables/columns\n")
    ## could be nice to say number of peaksData and chromData variables ?
  }
})

#' @importMethodsFrom S4Vectors [
#'
#' @importFrom MsCoreUtils i2index
#' @rdname hidden_aliases
setMethod("[", "ChromBackendMemory", function(x, i, j, ..., drop = FALSE) {
  i <- i2index(i, length = length(x))
  x@chromData <- x@chromData[i, ]
  x@peaksData <- x@peaksData[i]
  x
})

#' @rdname hidden_aliases
setMethod("$", "ChromBackendMemory", function(x, name) {
  if (name %in% union(chromVariables(x), names(coreChromVariables())))
    res <- chromData(x, columns = name, drop = TRUE)
  else if (name %in% peaksVariables(x))
    res <- peaksData(x, columns = name, drop = TRUE)
  else stop("The requested variable '", name, "' is not available")
  res
})

#' @rdname hidden_aliases
setReplaceMethod("$", "ChromBackendMemory", function(x, name, value) {
  if (length(x) && length(value) != length(x))
    stop("length of 'value' needs to match the number of chromatograms ",
         "in object.")
  if (name %in% peaksVariables(x)) {
    if (!is.list(value))
      stop("The value for peaksData should be a list")
    for (i in seq_along(value)) {
      x@peaksData[[i]][[name]] <- value[[i]]
      validPeaksData(x@peaksData)
    }
  } else {
    x@chromData[, name] <- value ## also this mean that new column will be created here, is that okay ?
    validChromData(x@chromData)
  }
  x
})
