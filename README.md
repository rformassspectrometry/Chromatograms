# Low level infrastructure to handle chromatographic data

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Chromatograms/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Chromatograms/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/Chromatograms/branch/main/graph/badge.svg?token=jy0Mid9gKn)](https://codecov.io/gh/rformassspectrometry/Chromatograms)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)


This package, part of the [*R for Mass
Spectrometry*](https://www.rformassspectrometry.org/) initiative,
provides support for handling and analysing chromatographic mass spectrometry
(MS) data.

It follows the API that is currenlty available in the
[`Chromatogram`](http://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
class in the [`MSnbase`](http://lgatto.github.io/MSnbase/index.html) package,
with support for different *backends* to provide/represent the data, similar to
the [`Spectra`](https://rformassspectrometry.github.io/Spectra/) package.


## General concept

A `Chromatograms` object is designed to contain multiple chromatographic data
(i.e. chromatogram entities). The data will be stored linearly, i.e. as a long
list of chromatograms. The `Chromatograms` object will be the main object for
the end user, providing functionality to access, filter or process
chromatographic data, with the actual chromatographic MS data being stored
within *backend* classes. Different implementations of backend classes can be
designed for high performance or low memory footprint.

The existing backend classes are: 

- `ChromBackendMemory`: a memory-based backend, storing the data in
  memory. This is the default backend and is used for testing purposes.
- The `ChromBackendMzR` inherits all slots and methods from the base
  `ChromBackendMemory` backend, providing additional functionality for reading 
  chromatographic data from mzML files.

- `ChromBackendSpectra`: The `ChromBackendSpectra` inherits all slots and
   methods from the base `ChromBackendMemory` backend, providing additional
   functionality for reading chromatographic data from `Spectra` objects.
   
These backend are then handled on a user level by the `Chromatograms` class,
which provides a unified interface to access and manipulate the chromatographic
data. 

## Installation

```
install.packages("BiocManager")
BiocManager::install("RforMassSpectrometry/Chromatograms")
```

## Contributions

Please the *R for Mass Spectrometry* [code of conduct](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct) and [contribution guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
