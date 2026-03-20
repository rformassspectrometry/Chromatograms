# Version 1.1

## Changes in 1.1.5

- Add `peakBoundary()` method for `Chromatograms` objects. Determines the
  retention time boundaries of the tallest peak in each chromatogram using
  `MsCoreUtils::valleys()` to locate flanking valleys, with a threshold-based
  fallback. Returns a matrix with `left_boundary` and `right_boundary` columns.

## Changes in 1.1.4

- `setBackend()` now clears the processing queue after switching backend,
  preventing queued processing steps from being applied twice (once during
  the data transfer and again on subsequent `peaksData()` calls).

- Fix `setBackend()` parallel branch (used for `ChromBackendMzR`) to correctly
  apply queued processing steps to each chunk before transferring data to
  the new backend.

- Major performance improvement in `.process_peaks_data()` for
  `ChromBackendSpectra`. Key optimizations include: pre-extracting `mz()` and
  `intensity()` as plain R lists (avoiding slow `SimpleNumericList` indexing),
  global retention time pre-filtering, a fast path for TIC/BPC cases, and
  using `findInterval()` with `cumsum()` for m/z range lookups. Combined,
 `setBackend()` showed 9x speed up for 1000 chromatograms.

- Add optimized `intensity()` and `rtime()` accessors for `ChromBackendMemory`
  using direct `[[` extraction instead of the slower `[, col, drop]` path
  through `peaksData()`.

- Further accessor optimizations: `intensity()`, `rtime()`, and `lengths()` on
  `Chromatograms` now bypass `peaksData()` dispatch when the processing queue
  is empty. `peaksData()` on `ChromBackendMemory` uses a fast `[[` path for
  single-column requests. Direct `lengths()` methods added for all backends
  using `nrow()` instead of going through `intensity()`.

- Replace `do.call(rbind, ...)` with `data.table::rbindlist()` in
  `chromExtract()` for `ChromBackendMemory`, `ChromBackendMzR`, and
  `ChromBackendSpectra` for faster row-binding of many data.frames.

- Replace `replicate(n, .EMPTY_PEAKS_DATA, simplify = FALSE)` with
  `rep(list(.EMPTY_PEAKS_DATA), n)` across backends to avoid repeated
  expression evaluation overhead.

## Changes in 1.1.3

- Add `filterEmptyChromatograms()` function to remove empty chromatograms
  (i.e., chromatograms without peaks) from a `Chromatograms` or
  `ChromBackend` object.

- Add `concatenateChromatograms()` function and `c()` method to combine
  multiple `Chromatograms` objects into a single object. Also add `split()`
  method to split a `Chromatograms` object based on a grouping factor.

- Add `extrapolate` parameter to `imputePeaksData()` (default `FALSE`).
  When `TRUE`, leading/trailing `NA` values outside the range of observed
  data are extrapolated. When `FALSE` (default), only interpolation is
  performed and edge `NA` values remain as `NA`.

## Changes in 1.1.2

- Fix `peaksData()` for `ChromBackendSpectra` to return data in the correct
  row order when multiple chromatograms share the same `chromSpectraIndex`.
  This bug caused `setBackend()` to produce mismatched `chromData` and
  `peaksData` when converting from `ChromBackendSpectra` to
  `ChromBackendMemory` with objects containing multiple EICs.

## Changes in 1.1.1

- Aligned the package with the Bioconductor 3.22 release.
- Expanded the vignette to cover ChromBackendSpectra usage, chromatogram
  extraction with `chromExtract()`, and imputation workflows via
  `imputePeaksData()`.
- Added `spectraSortIndex()` for `ChromBackendSpectra` to compute the desired
  retention-time order on demand, avoiding the need to keep on-disk `Spectra`
  objects sorted in memory.

# Version 0.99

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
