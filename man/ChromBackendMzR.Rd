% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChromBackendMzR.R
\name{ChromBackendMzR}
\alias{ChromBackendMzR}
\alias{backendInitialize,ChromBackendMzR-method}
\title{Chromatographic Data Backend for Reading mzML Files}
\usage{
ChromBackendMzR()

\S4method{backendInitialize}{ChromBackendMzR}(object, files = character(), BPPARAM = bpparam(), ...)
}
\arguments{
\item{object}{A \code{ChromBackendMzR} object.}

\item{files}{A character vector of file paths to mzML files.}

\item{BPPARAM}{Parallel setup configuration. See \code{\link[BiocParallel:register]{BiocParallel::bpparam()}}
for more information.}

\item{...}{Additional parameters to be passed.}
}
\description{
The \code{ChromBackendMzR} inherits all slots and methods from the base
\code{ChromBackendMemory} backend, providing additional functionality for reading
chromatographic data from mzML files.

Unlike the \code{ChromBackendMemory} backend, the \code{ChromBackendMzR} backend
should have the \emph{dataOrigin} chromatographic variables populated with the
file path of the mzML file from which the chromatographic data was read.

Note that the \code{ChromBackendMzR} backend is read-only and does not support
direct modification of chromatographic data. However, it does support
\code{peaksData} slot replacement, which will modify the \code{peaksData} slot but not
the local mzML files. This is indicated by the "inMemory" slot being set to
TRUE.

Implementing functionalities with the \code{ChromBackendMzR} backend should be
simplified as much as possible and reuse the methods already implemented for
\code{ChromBackendMemory} when possible.
}
\examples{
library(mzR)
library(msdata)

## Load an mzML file
MRM_file <- system.file("proteomics", "MRM-standmix-5.mzML.gz",
    package = "msdata"
)

## Initialize the ChromBackendMzR object
be_empty <- ChromBackendMzR()
be <- backendInitialize(be_empty, files = MRM_file, BPPARAM = SerialParam())

}
\author{
Philippine Louail
}
