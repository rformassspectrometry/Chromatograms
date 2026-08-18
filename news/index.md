# Changelog

## Version 1.3

### Changes in 1.3.3

- Major `ChromBackendSpectra` performance improvements:
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md),
  `Chromatograms(spectra)`,
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  and `[` no longer re-validate the wrapped `Spectra` on every call
  (which re-stated every backing file), so they no longer scale with the
  number of files.

- Improve performance of
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  for `ChromBackendSpectra` with overlapping chromatogram windows
  (e.g. from
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)):
  each spectrum is aggregated once and shared across the windows it
  falls in, instead of once per window. Results are unchanged; the
  speed-up grows with the number of overlapping windows.

- Improve performance of
  [`peakBoundary()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  valleys flanking the apex are located by scanning outwards from it
  rather than scanning the whole chromatogram. Results are unchanged.

- Order `dataOrigin` by first appearance when computing the spectra sort
  index, consistent with
  [`backendParallelFactor()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html).

### Changes in 1.3.2

- Change
  [`plotChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/plotChromatograms.md)
  and
  [`plotChromatogramsOverlay()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
  to methods.
- Import
  [`compareChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  from *ProtGenerics*.

### Changes in 1.3.1

- Addition of logo and align to BioC 3.23 release.

- Fix `ChromBackendSpectra` `spectraSortIndex` test to shuffle spectra
  first, avoiding a spurious failure when input data is already sorted.

## Version 1.1

### Changes in 1.1.8

- Improve performance of
  [`matchRtime()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md).

- Fix
  [`compareChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md):
  `...` arguments (e.g. `tolerance`) are now routed to `MAPFUN` or `FUN`
  based on their formal parameters, preventing errors when `FUN = cor`
  received unknown arguments.

### Changes in 1.1.7

- Improve performance of `.prepare_spectra_input()`: spectra are now
  pre-filtered to the non-overlapping union of EIC retention time ranges
  using
  [`MsCoreUtils::reduce()`](https://rdrr.io/pkg/MsCoreUtils/man/reduce.html)
  and
  [`Spectra::filterRanges()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  and peak data is loaded in a single
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  call instead of separate
  [`mz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  and
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  calls. This reduces I/O and memory usage, especially for file-backed
  backends.

### Changes in 1.1.6

- Add
  [`compareChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  and
  [`matchRtime()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  for pairwise similarity of chromatographic intensity profiles.

### Changes in 1.1.5

- Add
  [`peakBoundary()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  method for `Chromatograms` objects. Determines the retention time
  boundaries of the tallest peak in each chromatogram using
  [`MsCoreUtils::valleys()`](https://rdrr.io/pkg/MsCoreUtils/man/valleys.html)
  to locate flanking valleys, with a threshold-based fallback. Returns a
  matrix with `left_boundary` and `right_boundary` columns.

### Changes in 1.1.4

- [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  now clears the processing queue after switching backend, preventing
  queued processing steps from being applied twice (once during the data
  transfer and again on subsequent
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  calls).

- Fix
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  parallel branch (used for `ChromBackendMzR`) to correctly apply queued
  processing steps to each chunk before transferring data to the new
  backend.

- Major performance improvement in `.process_peaks_data()` for
  `ChromBackendSpectra`. Key optimizations include: pre-extracting
  [`mz()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  and
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  as plain R lists (avoiding slow `SimpleNumericList` indexing), global
  retention time pre-filtering, a fast path for TIC/BPC cases, and using
  [`findInterval()`](https://rdrr.io/r/base/findInterval.html) with
  [`cumsum()`](https://rdrr.io/r/base/cumsum.html) for m/z range
  lookups. Combined,
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  showed 9x speed up for 1000 chromatograms.

- Add optimized
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  and
  [`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  accessors for `ChromBackendMemory` using direct `[[` extraction
  instead of the slower `[, col, drop]` path through
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md).

- Further accessor optimizations:
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  [`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  and [`lengths()`](https://rdrr.io/r/base/lengths.html) on
  `Chromatograms` now bypass
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  dispatch when the processing queue is empty.
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  on `ChromBackendMemory` uses a fast `[[` path for single-column
  requests. Direct [`lengths()`](https://rdrr.io/r/base/lengths.html)
  methods added for all backends using
  [`nrow()`](https://rdrr.io/pkg/BiocGenerics/man/nrow.html) instead of
  going through
  [`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

- Replace `do.call(rbind, ...)` with
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
  in
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  for `ChromBackendMemory`, `ChromBackendMzR`, and `ChromBackendSpectra`
  for faster row-binding of many data.frames.

- Replace `replicate(n, .EMPTY_PEAKS_DATA, simplify = FALSE)` with
  `rep(list(.EMPTY_PEAKS_DATA), n)` across backends to avoid repeated
  expression evaluation overhead.

### Changes in 1.1.3

- Add
  [`filterEmptyChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  function to remove empty chromatograms (i.e., chromatograms without
  peaks) from a `Chromatograms` or `ChromBackend` object.

- Add
  [`concatenateChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/concatenateChromatograms.md)
  function and [`c()`](https://rdrr.io/r/base/c.html) method to combine
  multiple `Chromatograms` objects into a single object. Also add
  [`split()`](https://rdrr.io/r/base/split.html) method to split a
  `Chromatograms` object based on a grouping factor.

- Add `extrapolate` parameter to
  [`imputePeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  (default `FALSE`). When `TRUE`, leading/trailing `NA` values outside
  the range of observed data are extrapolated. When `FALSE` (default),
  only interpolation is performed and edge `NA` values remain as `NA`.

### Changes in 1.1.2

- Fix
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  for `ChromBackendSpectra` to return data in the correct row order when
  multiple chromatograms share the same `chromSpectraIndex`. This bug
  caused
  [`setBackend()`](https://rdrr.io/pkg/ProtGenerics/man/backendInitialize.html)
  to produce mismatched `chromData` and `peaksData` when converting from
  `ChromBackendSpectra` to `ChromBackendMemory` with objects containing
  multiple EICs.

### Changes in 1.1.1

- Aligned the package with the Bioconductor 3.22 release.
- Expanded the vignette to cover ChromBackendSpectra usage, chromatogram
  extraction with
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md),
  and imputation workflows via
  [`imputePeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md).
- Added `spectraSortIndex()` for `ChromBackendSpectra` to compute the
  desired retention-time order on demand, avoiding the need to keep
  on-disk `Spectra` objects sorted in memory.

## Version 0.99

### Changes in 0.99.7

- Add
  [`chromExtract()`](https://rformassspectrometry.github.io/Chromatograms/reference/Chromatograms.md)
  method to generate a new `Chromatograms` object from an existing one
  by extracting a subset of chromatograms based on retention times
  (optionally m/z) boundaries.
- Add
  [`imputePeaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  method to impute missing values in the chromatographic peaks data.
- Fix
  [`factorize()`](https://rformassspectrometry.github.io/Chromatograms/reference/ChromBackend.md)
  so that the parameter `factorize.by` can take a `character` vector of
  length 1.

### Changes in 0.99.6

- Add `Spectra` dependency.

### Changes in 0.99.5

- Add `IRanges` dependency

### Changes in 0.99.4

- Add dependencies to Vignette.

### Changes in 0.99.0

- General documentation and formatting fixes for Bioconductor
  submission.

### Changes in 0.6.0

- Addition of `ChromBackendSpectra` class and its respective methods.
- Addition of
  [`plotChromatograms()`](https://rformassspectrometry.github.io/Chromatograms/reference/plotChromatograms.md)
  and
  [`plotChromatogramsOverlay()`](https://rformassspectrometry.github.io/Chromatograms/reference/hidden_aliases.md)
  functions.
- Addition of the `extractByIndex` implementation in the backends.

### Changes in 0.5.0

- Addition of `ChromBackendMzR` and its respective methods.
- Addition of the Chromatograms vignette, which provides an overview of
  the object and related functionalities.

### Changes in 0.4.0

- Addition of
  [`peaksData()`](https://rformassspectrometry.github.io/Chromatograms/reference/peaksData.md)
  and implementation of chunkwise (and therefore paralleled) processing
  of `Chromatograms` object.
- Addition of
  [`addProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html),
  [`applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html),
  [`processingChunkFactor()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html),
  and
  [`processingChunkSize()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html).

### Changes in 0.3.0

- Addition of
  [`filterChromData()`](https://rformassspectrometry.github.io/Chromatograms/reference/chromData.md)
  method for `ChromBackend`.
- Creation of the `Chromatograms` class and implementation of basic
  accessor methods.
- Addition of basic plotting functions.

### Changes in 0.2.0

- Addition of `ChomBackendMemory` class and associated methods

### Changes in 0.1.0

- Addition of basic `ChromBackend` class and default methods
