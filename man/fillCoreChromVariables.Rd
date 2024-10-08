% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChromBackend-functions.R
\name{fillCoreChromVariables}
\alias{fillCoreChromVariables}
\title{Fill data frame with columns for missing core chrom variables}
\usage{
fillCoreChromVariables(x = data.frame())
}
\arguments{
\item{x}{\code{data.frame} with potentially present core chrom variable columns}
}
\value{
input data frame \code{x} with missing core variables added (with the
correct data type).
}
\description{
\code{fillCoreChromVariables()} fills a provided \code{data.frame}
with columns for eventually missing \emph{core} chromatogram variables.
The missing core variables are added as new columns with missing values
(\code{NA}) of the correct data type.
Use \code{\link[=coreChromVariables]{coreChromVariables()}} to list the set of core variables and their data
types.
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
