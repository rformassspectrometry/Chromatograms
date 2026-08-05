# Low level infrastructure to handle chromatographic data

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![years in
bioc](http://bioconductor.org/shields/years-in-bioc/Chromatograms.svg)](https://bioconductor.org/packages/release/bioc/html/Chromatograms.html)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/Chromatograms/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/Chromatograms/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/Chromatograms/branch/main/graph/badge.svg?token=jy0Mid9gKn)](https://codecov.io/gh/rformassspectrometry/Chromatograms)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

This package, part of the [*R for Mass
Spectrometry*](https://www.rformassspectrometry.org/) initiative,
provides a powerful and expandable infrastructure for handling and
analysing chromatographic mass spectrometry (MS) data.

It will replace the legacy classes to handle chromatographic data in R
provided by the [*MSnbase*](http://lgatto.github.io/MSnbase/index.md)
package.

This package is part of **Bioconductor**:
<https://bioconductor.org/packages/Chromatograms>.

## ⚙️ General concept

A `Chromatograms` object is designed to contain multiple chromatographic
data (i.e. chromatogram entities). The data will be stored linearly,
i.e. as a long list of chromatograms. The `Chromatograms` object will be
the main object for the end user, providing functionality to access,
filter or process chromatographic data, with the actual chromatographic
MS data being stored within *backend* classes. Different implementations
of backend classes can be designed for high performance or low memory
footprint.

The existing backend classes are:

- `ChromBackendMemory`: a memory-based backend, storing the data in
  memory. This is the default backend and is used for testing purposes.

- The `ChromBackendMzR` inherits all slots and methods from the base
  `ChromBackendMemory` backend, providing additional functionality for
  reading chromatographic data from mzML files.

- `ChromBackendSpectra`: The `ChromBackendSpectra` inherits all slots
  and methods from the base `ChromBackendMemory` backend, providing
  additional functionality for reading chromatographic data from
  `Spectra` objects.

These backend are then handled on a user level by the `Chromatograms`
class, which provides a unified interface to access and manipulate the
chromatographic data.

## ⤵️ Installation

    install.packages("BiocManager")
    BiocManager::install("Chromatograms")

## 🤝 Contribution

Please help us improving and completing the package! Any type of
contribution welcome 👐 - including discussions, suggestions or actual
code. Don’t be afraid - we’re friendly ☺️! 👉 get involved by opening an
[issue](https://github.com/rformassspectrometry/Chromatograms/issues).

Please also check out the [**RforMassSpectrometry Contributions
Guide**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).

### 📜 Code of Conduct

We follow the [**RforMassSpectrometry Code of
Conduct**](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#code-of-conduct)
to maintain an inclusive and respectful community.

## License

This package is licensed under the **Artistic 2.0** license: 📄
<https://opensource.org/license/Artistic-2.0>

Documentation (manuals, vignettes) is licensed under **CC BY-NC-SA
4.0**: 📄 <https://creativecommons.org/licenses/by-nc-sa/4.0/>

------------------------------------------------------------------------

# Funding information

Part of this work was funded by the **European Union** under the
**HORIZON-MSCA-2021** project **101073062: HUMAN – Harmonising and
Unifying Blood Metabolic Analysis Networks**.

![EU
Logo](https://github.com/rformassspectrometry/Metabonaut/raw/main/vignettes/images/EULogo.jpg)

EU Logo
