library("testthat")
library("Chromatograms")

# A data.frame with chromatogram variables.
cdata <- data.frame(msLevel = c(1L, 1L, 1L),
                    mz = c(112.2, 123.3, 134.4),
                    chromIndex = c(1L, 2L, 3L))

# Retention time and intensity values for each chromatogram.
pdata <- list(
    data.frame(rtime = c(12.4, 12.8, 13.2, 14.6),
               intensity = c(123.3, 153.6, 2354.3, 243.4)),
    data.frame(rtime = c(45.1, 46.2),
               intensity = c(100, 80.1)),
    data.frame(rtime = c(12.4, 12.8, 13.2, 14.6),
               intensity = c(123.3, 153.6, 2354.3, 243.4))
)

be_empty <- new("ChromBackendMemory")
be_cd <- backendInitialize(be_empty, chromData = cdata)
be <- backendInitialize(be_empty, chromData = cdata, peaksData = pdata)

## Run tests with the unit test suite defined in the Chromatograms package to
## ensure compliance with the definitions of the ChromBackend interface/class.
test_suite <- system.file("test_backends", "test_ChromBackend",
                          package = "Chromatograms")
test_dir(test_suite, stop_on_failure = TRUE)

c_empty <- Chromatograms()
c_full <- Chromatograms(be)

test_check("Chromatograms")
