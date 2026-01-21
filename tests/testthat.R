library(testthat)
library(Chromatograms)
library(Spectra)
library(MsBackendMetaboLights) ## remove MsbackendMetaboLights dependency for test and examples
library(IRanges)
library(S4Vectors)

### Test ChromBackendSpectra
be <- backendInitialize(MsBackendMetaboLights(),
    mtblsId = "MTBLS39",
    filePattern = c("63B.cdf")
)
s <- Spectra(be)
s <- setBackend(s, MsBackendMemory())
be_empty <- new("ChromBackendSpectra")
be <- backendInitialize(be_empty, s)
test_suite <- system.file("test_backends", "test_ChromBackend",
    package = "Chromatograms"
)
test_dir(test_suite, stop_on_failure = TRUE)
be_sp <- be
c_sp <- Chromatograms(be)

### Test ChrombackendMzR
# fetch files
MRM_file <- system.file("proteomics", "MRM-standmix-5.mzML.gz",
    package = "msdata"
)
be_empty <- ChromBackendMzR()
be <- backendInitialize(be_empty, files = MRM_file, BPPARAM = SerialParam())

test_suite <- system.file("test_backends", "test_ChromBackend",
    package = "Chromatograms"
)
test_dir(test_suite, stop_on_failure = TRUE)

be_mzr <- be
c_mzr <- Chromatograms(be)

### Test ChrombackendMemory
cdata <- data.frame(
    msLevel = c(1L, 1L, 1L),
    mz = c(112.2, 123.3, 134.4),
    dataOrigin = c("mem1", "mem2", "mem3")
)
pdata <- list(
    data.frame(
        rtime = c(12.4, 12.8, 13.2, 14.6),
        intensity = c(123.3, 153.6, 2354.3, 243.4)
    ),
    data.frame(
        rtime = c(45.1, 46.2),
        intensity = c(100, 80.1)
    ),
    data.frame(
        rtime = c(12.4, 12.8, 13.2, 14.6),
        intensity = c(123.3, 153.6, 2354.3, 243.4)
    )
)

be_empty <- new("ChromBackendMemory")
be_cd <- backendInitialize(be_empty, chromData = cdata)
be <- backendInitialize(be_empty, chromData = cdata, peaksData = pdata)
test_dir(test_suite, stop_on_failure = TRUE)

c_empty <- Chromatograms()
c_full <- Chromatograms(be)

test_check("Chromatograms")
