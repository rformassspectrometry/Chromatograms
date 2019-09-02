library("testthat")
library("Chromatograms")

fl <- msdata::proteomics(full.names = TRUE)[1]
mrm_mzr <- backendInitialize(ChromBackendMzR(), fl)

test_check("Chromatograms")
