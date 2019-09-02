#' @include hidden_aliases.R
NULL

.valid_chrom_backend_files_exist <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) && !all(file.exists(x)))
        return(paste0("File(s) ", paste(x[!file.exists(x)], collapse = ", "),
                      " not found"))
    NULL
}

#' @rdname ChromBackend
#'
#' @export ChromBackendMzR
ChromBackendMzR <- function() {
    if (!requireNamespace("mzR", quietly = TRUE))
        stop("The use of 'ChromBackendMzR' requires package 'mzR'. Please ",
             "install with 'Biobase::install(\"mzR\")'")
    new("ChromBackendMzR")
}

#' Read the chromatogram header for each chromatogram from the MS file `x`
#'
#' @author Johannes Rainer
#'
#' @return `DataFrame` with the header.
#'
#' @noRd
.mzR_chrom_header <- function(x = character()) {
    if (length(x) != 1)
        stop("'x' should have length 1")
    requireNamespace("mzR", quietly = TRUE)
    msd <- mzR::openMSfile(x)
    on.exit(mzR::close(msd))
    hdr <- mzR::chromatogramHeader(msd)
    colnames(hdr)[colnames(hdr) == "chromatogramIndex"] <- "chromIndex"
    colnames(hdr)[colnames(hdr) == "precursorCollisionEnergy"] <- "collisionEnergy"
    colnames(hdr)[colnames(hdr) == "precursorIsolationWindowTargetMZ"] <- "precursorMz"
    colnames(hdr)[colnames(hdr) == "precursorIsolationWindowLowerOffset"] <- "precursorMzMin"
    colnames(hdr)[colnames(hdr) == "precursorIsolationWindowUpperOffset"] <- "precursorMzMax"
    colnames(hdr)[colnames(hdr) == "productIsolationWindowTargetMZ"] <- "productMz"
    colnames(hdr)[colnames(hdr) == "productIsolationWindowLowerOffset"] <- "productMzMin"
    colnames(hdr)[colnames(hdr) == "productIsolationWindowUpperOffset"] <- "productMzMax"
    hdr$precursorMzMin <- hdr$precursorMz - hdr$precursorMzMin
    hdr$precursorMzMax <- hdr$precursorMz + hdr$precursorMzMax
    hdr$productMzMin <- hdr$productMz - hdr$productMzMin
    hdr$productMzMin <- hdr$productMz + hdr$productMzMax
    DataFrame(hdr)
}

#' Read chromatogram data from a single mzML file.
#'
#' @param x `character(1)` with the file to read from.
#'
#' @param chromIndex (required) indices of chromatograms from which the data
#'     should be retrieved.
#'
#' @author Johannes Rainer
#'
#' @return `list` of `data.frame` with columns `"rtime"` and `"intensity"`
#'
#' @noRd
.mzR_chromatograms <- function(x = character(), chromIndex = integer()) {
    if (length(x) != 1)
        stop("'x' should have length 1")
    msd <- mzR::openMSfile(x)
    on.exit(mzR::close(msd))
    chrs <- mzR::chromatogram(msd, chromIndex)
    if (is.data.frame(chrs))
        chrs <- list(chrs)
    lapply(chrs, function(z) {
        colnames(z) <- c("rtime", "intensity")
        as.matrix(z)
    })
}

#' @importFrom IRanges NumericList
#'
#' @description
#'
#' Helper to build the chromData `DataFrame` for `ChromBackendMzR` backend.
#'
#' @param x `ChromBackend` with a `chromData` slot and the `pairs` method
#'     defined.
#'
#' @param columns `character` defining the columns to return.
#'
#' @noRd
.chrom_data_mzR <- function(x, columns) {
    cn <- colnames(x@chromData)
    if(!nrow(x@chromData)) {
        res <- lapply(
            .CHROMATOGRAMS_DATA_COLUMNS[!(names(.CHROMATOGRAMS_DATA_COLUMNS)
                %in% c("rtime", "intensity"))], do.call, args = list())
        res <- DataFrame(res)
        res$rtime <- NumericList(compress = FALSE)
        res$intensity <- NumericList(compress = FALSE)
        return(res[, columns, drop = FALSE])
    }
    not_found <- setdiff(columns, c(cn, names(.CHROMATOGRAMS_DATA_COLUMNS)))
    if (length(not_found))
        stop("Column(s) ", paste(not_found, collapse = ", "),
             " not available")
    cols <- columns[columns %in% cn]
    res <- asVectorDataFrame(x@chromData[, cols, drop = FALSE])
    any_rtime <- any(columns == "rtime")
    any_int <- any(columns == "intensity")
    if (any_rtime || any_int) {
        prs <- pairs(x)
        if (any_rtime)
            res$rtime <- NumericList(lapply(prs, "[", , 1), compress = FALSE)
        if (any_int)
            res$intensity <- NumericList(lapply(prs, "[", , 2),
                                         compress = FALSE)
    }
    other_cols <- setdiff(columns, c(cols, "rtime", "intensity"))
    if (length(other_cols)) {
        other_res <- lapply(other_cols, .get_rle_column, x = x@chromData)
        names(other_res) <- other_cols
        res <- cbind(res, as(other_res, "DataFrame"))
    }
    res[, columns, drop = FALSE]
}
