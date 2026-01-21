# Version 1.1.1

- Aligned the package with the Bioconductor 3.22 release.
- Expanded the vignette to cover ChromBackendSpectra usage, chromatogram
  extraction with `chromExtract()`, and imputation workflows via
  `imputePeaksData()`.
- Added `spectraSortIndex()` for `ChromBackendSpectra` to compute the desired
  retention-time order on demand, avoiding the need to keep on-disk `Spectra`
  objects sorted in memory.

# Version 0.99.7

## Changes in 0.99.7

- Add `chromExtract()` method to generate a new `Chromatograms` object from an
  existing one by extracting a subset of chromatograms based on retention 
  times (optionally m/z) boundaries.
- Add `imputePeaksData()` method to impute missing values in the 
  chromatographic peaks data.
- Fix `factorize()` so that the parameter `factorize.by` can take a `character`
  vector of length 1. 

## Changes in 0.99.6

- Add `Spectra` dependency.

## Changes in 0.99.5

- Add `IRanges` dependency

## Changes in 0.99.4

- Add dependencies to Vignette.

## Changes in 0.99.0

- General documentation and formatting fixes for Bioconductor submission.

## Changes in 0.6.0

- Addition of `ChromBackendSpectra` class and its respective methods.
- Addition of `plotChromatograms()` and `plotChromatogramsOverlay()` functions.
- Addition of the `extractByIndex` implementation in the backends.

## Changes in 0.5.0

- Addition of `ChromBackendMzR` and its respective methods.
- Addition of the Chromatograms vignette, which provides an overview of the
  object and related functionalities.

## Changes in 0.4.0
- Addition of `peaksData()` and implementation of chunkwise (and therefore
  paralleled) processing of `Chromatograms` object.
- Addition of `addProcessing()`, `applyProcessing()`, `processingChunkFactor()`,
  and `processingChunkSize()`.

## Changes in 0.3.0
- Addition of `filterChromData()` method for `ChromBackend`.
- Creation of the `Chromatograms` class and implementation of basic accessor
  methods.
- Addition of basic plotting functions.

## Changes in 0.2.0
- Addition of `ChomBackendMemory` class and associated methods

## Changes in 0.1.0
- Addition of basic `ChromBackend` class and default methods
