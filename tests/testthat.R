library("testthat")
library("Chromatograms")

fl <- msdata::proteomics(full.names = TRUE)[1]

test_check("Chromatograms")
