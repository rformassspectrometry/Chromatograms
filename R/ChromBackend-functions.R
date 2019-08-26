#' @include hidden_aliases.R
NULL

.valid_chrom_backend_data_storage <- function(x) {
    if (anyNA(x))
        return("'NA' values in dataStorage are not allowed.")
    NULL
}
