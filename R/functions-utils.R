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
