% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Chromatograms.R
\name{Chromatograms}
\alias{Chromatograms}
\alias{Chromatograms-class}
\alias{[,Chromatograms-method}
\alias{[<-,Chromatograms-method}
\alias{[[,Chromatograms-method}
\alias{[[<-,Chromatograms-method}
\alias{Chromatograms,ChromBackendOrMissing-method}
\alias{Chromatograms,Spectra-method}
\alias{setBackend,Chromatograms,ChromBackend-method}
\alias{$,Chromatograms-method}
\alias{$<-,Chromatograms-method}
\alias{factorize,Chromatograms-method}
\title{The Chromatograms class to manage and access chromatographic data}
\usage{
\S4method{Chromatograms}{ChromBackendOrMissing}(object = ChromBackendMemory(), processingQueue = list(), ...)

\S4method{Chromatograms}{Spectra}(
  object,
  summarize.method = c("sum", "max"),
  chromData = data.frame(),
  factorize.by = c("msLevel", "dataOrigin"),
  ...
)

\S4method{setBackend}{Chromatograms,ChromBackend}(
  object,
  backend,
  f = processingChunkFactor(object),
  BPPARAM = SerialParam(),
  ...
)

\S4method{$}{Chromatograms}(x, name)

\S4method{$}{Chromatograms}(x, name) <- value

\S4method{[}{Chromatograms}(x, i, j, ..., drop = FALSE)

\S4method{[[}{Chromatograms}(x, i, j, ...)

\S4method{[[}{Chromatograms}(x, i, j, ...) <- value

\S4method{factorize}{Chromatograms}(object, factorize.by = c("msLevel", "dataOrigin"), ...)
}
\arguments{
\item{object}{A \link{Chromatograms} object.}

\item{processingQueue}{\link{list} a list of processing steps (i.e. functions) to
be applied to the chromatographic data. The processing steps are
applied in the order they are listed in the \code{processingQueue}.}

\item{...}{Additional arguments.}

\item{summarize.method}{For Chromatograms created with a \code{Spectra} object:
A \code{character} vector with the name of the function to be used to
summaries the spectra data intensity. The available methods are "sum"
and "max". The default is "sum".}

\item{chromData}{For \code{Chromatograms()} build from a \code{Spectra} object backend,
a \code{data.frame} with the chromatographic data. If not provided
(or if empty), a default \code{data.frame} with the core chromatographic
variables will be created.}

\item{factorize.by}{A \code{character} vector with the names of the variables in
the \code{Spectra} object and the \code{chromData} slot that should be used
to factorize the \code{Spectra} object data to generate the
chromatographic data.}

\item{backend}{\link{ChromBackend} object providing the raw data for the
\code{Chromatograms} object.}

\item{f}{\code{factor} defining the grouping to split the \code{Chromatograms} object.}

\item{BPPARAM}{Parallel setup configuration. See \code{\link[BiocParallel:register]{BiocParallel::bpparam()}}
for more information.}

\item{x}{A \link{Chromatograms} object.}

\item{name}{A \code{character} string specifying the name of the variable to
access.}

\item{value}{The value to replace the variable with.}

\item{i}{For \code{[}: \code{integer}, \code{logical} or \code{character} to subset the object.}

\item{j}{For \code{[} and \code{[[}: ignored.}

\item{drop}{For \code{[}: \code{logical(1)} default to \code{FALSE}.}
}
\description{
The \code{Chromatograms} class encapsules chromatographic data and related
metadata. The chromatographic data is represented by a \emph{backend} extending
the virtual \link{ChromBackend} class which provides the raw data to the
\code{Chromatograms} object. Different backends and their properties are
decribed in the \link{ChromBackend} class documentation.
}
\note{
This needs to be discussed, if we want for example to be able to set a
a backend to \code{ChromBackendMzR} we need to implement backendInitialize()
better. = Support peaksData and chromData as arguments AND have a way to
write .mzml files (which we do not have for chromatographic data).
}
\section{Creation of objects}{


\code{Chromatograms} objects can be created using the \code{Chromatograms()}
construction function. Either by providing a \code{ChromBackend} object or by
providing a \code{Spectra} object. The \code{Spectra} object will be used to generate
a \code{Chromatograms} object with a backend of class \code{\link{ChromBackendSpectra}}.
}

\section{Data stored in a \code{Chromatograms} object}{


The \code{Chromatograms} object is a container for chromatographic data, which
includes peaks data (\emph{retention time} and related intensity values, also
referred to as \emph{peaks data variables} in the context of \code{Chromatograms}) and
metadata of individual chromatogram (so called \emph{chromatograms variables}).
While a core set of chromatograms variables (the
\code{coreChromatogramsVariables()}) and peaks data variables (the
\code{corePeaksVariables()}) are guaranteed to be provided by a \code{Chromatograms},
it is possible to add arbitrary variables to a \code{Chromatograms} object.

The \code{Chromatograms} object is designed to contain chromatographic data of a
(large) set of chromatograms. The data is organized \emph{linearly} and can be
thought of a list of chromatograms, i.e. each element in the \code{Chromatograms}
is one chromatogram.

The \emph{chromatograms variables} information in the \code{Chromatograms} object can
be accessed using the \code{chromData()} function. Specific chromatograms
variables can be accessed by either precising the \code{"columns"} parameter in
\code{chromData()} or using \code{$}. \code{chromData} can be accessed, replaced but
also filtered/subsetted. Refer to the \link{chromData} documentation for more
details.

The \emph{peaks data variables} information in the \code{Chromatograms} object can be
accessed using the \code{peaksData()} function. Specific peaks variables can be
accessed by either precising the \code{"columns"} parameter in \code{peaksData()} or
using \code{$}. \code{peaksData} can be accessed, replaced but also
filtered/subsetted. Refer to the \link{peaksData} documentation for more details.
}

\section{Processing of \code{Chromatograms} objects}{


Functions that process the chromatograms data in some ways can be applied to
the object either directly or by using the \code{processingQueue} mechanism. The
\code{processingQueue} is a list of processing steps that are stored within the
object and only applied when needed. This was created so that the data can be
processed in a single step and is very useful for larger datasets. This is
even more true as this processing queue will call function that can be
applied on the data in a chunk-wise manner. This allows for parallel
processing of the data and reduces the memory demand. To read more about the
\code{processingQueue}, and how to parallelize your processes, see the
\link{processingQueue} documentation.
}

\examples{

library(MsBackendMetaboLights)
library(Spectra)

## Create a Chromatograms object from a Spectra object.
be <- backendInitialize(MsBackendMetaboLights(),
    mtblsId = "MTBLS39",
    filePattern = c("63B.cdf")
)
s <- Spectra(be)
be <- backendInitialize(new("ChromBackendSpectra"), s)
chr <- Chromatograms(be)

## Subset
chr[1:2]

## access a specific variables
chr[["msLevel"]]
chr$msLevel

## Replace data of a specific variable
chr$msLevel <- c(2L, 2L, 2L)

## Can re factorize the data
chr <- factorize(chr)

## Can also change the backend into memory
chr <- setBackend(chr, ChromBackendMemory())

chr

}
\seealso{
\link{chromData} for a general description of the chromatographic
metadata available in the object, as well as how to access, replace
and subset them.
\link{peaksData} for a general description of the chromatographic peaks
data available in the object, as well as how to access, replace and
subset them.
\link{processingQueue} for more information on the queuing of
processings and parallelization for larger dataset.
}
