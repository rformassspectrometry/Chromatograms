% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/ChromBackend-functions.R,
%   R/hidden_aliases.R, R/ChromBackendMemory.R, R/ChromBackendMzR.R,
%   R/ChromBackendSpectra.R, R/Chromatograms.R
\docType{methods}
\name{reset}
\alias{reset}
\alias{fillCoreChromVariables}
\alias{validChromData}
\alias{validPeaksData}
\alias{hidden_aliases}
\alias{backendMerge,ChromBackendMemory-method}
\alias{chromData,ChromBackendMemory-method}
\alias{chromData<-,ChromBackendMemory-method}
\alias{chromVariables,ChromBackendMemory-method}
\alias{peaksData,ChromBackendMemory-method}
\alias{peaksData<-,ChromBackendMemory-method}
\alias{peaksVariables,ChromBackendMemory-method}
\alias{isReadOnly,ChromBackendMemory-method}
\alias{show,ChromBackendMemory-method}
\alias{supportsSetBackend,ChromBackendMemory-method}
\alias{[,ChromBackendMemory-method}
\alias{$,ChromBackendMemory-method}
\alias{$<-,ChromBackendMemory-method}
\alias{show,ChromBackendMzR-method}
\alias{backendParallelFactor,ChromBackendMzR-method}
\alias{isReadOnly,ChromBackendMzR-method}
\alias{peaksData,ChromBackendMzR-method}
\alias{peaksData<-,ChromBackendMzR-method}
\alias{chromData<-,ChromBackendMzR-method}
\alias{supportsSetBackend,ChromBackendMzR-method}
\alias{[,ChromBackendMzR-method}
\alias{show,ChromBackendSpectra-method}
\alias{factorize,ChromBackendSpectra-method}
\alias{backendParallelFactor,ChromBackendSpectra-method}
\alias{isReadOnly,ChromBackendSpectra-method}
\alias{peaksData,ChromBackendSpectra-method}
\alias{peaksData<-,ChromBackendSpectra-method}
\alias{chromData<-,ChromBackendSpectra-method}
\alias{supportsSetBackend,ChromBackendSpectra-method}
\alias{[,ChromBackendSpectra-method}
\alias{show,Chromatograms-method}
\title{Fill data.frame with columns for missing core chromatogram variables.}
\usage{
reset(object, ...)

fillCoreChromVariables(x = data.frame())

validChromData(x = data.frame(), error = TRUE)

validPeaksData(x = list(), error = TRUE)

\S4method{backendMerge}{ChromBackendMemory}(object, ...)

\S4method{chromData}{ChromBackendMemory}(object, columns = chromVariables(object), drop = FALSE)

\S4method{chromData}{ChromBackendMemory}(object) <- value

\S4method{chromVariables}{ChromBackendMemory}(object)

\S4method{peaksData}{ChromBackendMemory}(object, columns = peaksVariables(object), drop = FALSE, ...)

\S4method{peaksData}{ChromBackendMemory}(object) <- value

\S4method{peaksVariables}{ChromBackendMemory}(object)

\S4method{isReadOnly}{ChromBackendMemory}(object)

\S4method{show}{ChromBackendMemory}(object)

\S4method{supportsSetBackend}{ChromBackendMemory}(object, ...)

\S4method{[}{ChromBackendMemory}(x, i, j, ..., drop = FALSE)

\S4method{$}{ChromBackendMemory}(x, name)

\S4method{$}{ChromBackendMemory}(x, name) <- value

\S4method{show}{ChromBackendMzR}(object)

\S4method{backendParallelFactor}{ChromBackendMzR}(object, ...)

\S4method{isReadOnly}{ChromBackendMzR}(object)

\S4method{peaksData}{ChromBackendMzR}(
  object,
  columns = peaksVariables(object),
  drop = FALSE,
  BPPARAM = SerialParam(),
  ...
)

\S4method{peaksData}{ChromBackendMzR}(object) <- value

\S4method{chromData}{ChromBackendMzR}(object) <- value

\S4method{supportsSetBackend}{ChromBackendMzR}(object, ...)

\S4method{[}{ChromBackendMzR}(x, i, j, ..., drop = TRUE)

\S4method{show}{ChromBackendSpectra}(object)

\S4method{factorize}{ChromBackendSpectra}(object, factorize.by = c("msLevel", "dataOrigin"), ...)

\S4method{backendParallelFactor}{ChromBackendSpectra}(object, ...)

\S4method{isReadOnly}{ChromBackendSpectra}(object)

\S4method{peaksData}{ChromBackendSpectra}(object, columns = peaksVariables(object), drop = FALSE, ...)

\S4method{peaksData}{ChromBackendSpectra}(object) <- value

\S4method{chromData}{ChromBackendSpectra}(object) <- value

\S4method{supportsSetBackend}{ChromBackendSpectra}(object, ...)

\S4method{[}{ChromBackendSpectra}(x, i, j, ..., drop = TRUE)

\S4method{show}{Chromatograms}(object)
}
\arguments{
\item{object}{A \link{Chromatograms} object.}

\item{x}{\code{list} representing the peaks data of a \code{Chromatograms}}

\item{error}{\code{logical(1)} whether an error should be thrown (the default)
if one or more columns don't have the correct data type.}
}
\value{
input data frame \code{x} with missing core variables added (with the
correct data type).

If the core variables have all the correct data type: an empty character.
If one or more core variables (columns) have the wrong data type the
function either throws an error (with \code{error = TRUE}) or returns a
\code{character} specifying which variables/columns don't have the correct
type (for \code{error = FALSE}).

Not applicable
}
\description{
\code{fillCoreChromVariables()} fills a provided \code{data.frame}
with columns for eventually missing \emph{core} chromatogram variables.
The missing core variables are added as new columns with missing values
(\code{NA}) of the correct data type.
Use \code{\link[=coreChromVariables]{coreChromVariables()}} to list the set of core variables and their data
types.

\code{validChromData()} checks that columns, representing \emph{core} chromatogram
variables are of the correct data type.

For S4 methods that require a documentation entry but only clutter the
index.

This method returns the chromatographic data stored in the
backend. If not specified otherwise it will return all defined columns in
the chromData slot as well as adding the \code{coreChromVariables} missing with
NA values.
}
\examples{

## Define a data frame
a <- data.frame(msLevel = c(1L, 1L, 2L), other_column = "b")

## Add missing core chromatogram variables to this data frame
fillCoreChromVariables(a)

## The data.frame thus contains columns for all core chromatogram
## variables in the respective expected data type (but filled with
## missing values).

}
\keyword{internal}
