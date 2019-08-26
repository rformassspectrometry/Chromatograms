#' Helper function to check the data type even if `x` is an `Rle`.
#'
#' @param x
#'
#' @param class2 `character(1)` specifying the class for which should be tested.
#'
#' @author Johannes Rainer
#'
#' @noRd
.is_class <- function(x, class2) {
    if (is(x, "Rle"))
        is(x@values, class2)
    else is(x, class2)
}

#' Helper to convert i to integer index for subsetting.
#'
#' @param i `character` `logical` or `integer` used in `[i]` for subsetting
#'
#' @param length `integer` representing the `length` of the object to be
#'     subsetted
#'
#' @param names `character` with the names (rownames or similar) of the object
#'
#' @return `integer` with the indices
#'
#' @author Johannes Rainer
#'
#' @noRd
.i_to_index <- function(i, length, names = NULL) {
    if (is.character(i)) {
        if (!length(names))
            stop("can not subset by name, object does not have names")
        if (!all(i %in% names))
            stop("not all names present in object")
        i <- match(i, names)
    }
    if (is.logical(i)) {
        if (length(i) != length)
            stop("if i is logical it has to match the length of the object (",
                 length, ")")
        i <- which(i)
    }
    if (is.numeric(i))
        i <- as.integer(i)
    if (is.integer(i) && !all(abs(i) %in% seq_len(length)))
        stop("index out of bounds: index has to be between 1 and ", length)
    i
}
