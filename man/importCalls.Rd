% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{importCalls}
\alias{importCalls}
\title{Load structural variant calls}
\usage{
importCalls(file = NULL, ID = NULL, remove.uncertain = TRUE)
}
\arguments{
\item{file}{A file containing structural variant calls.}

\item{ID}{A metadata column to be appended.}

\item{remove.uncertain}{Set to \code{TRUE} if calls with question mark should be removed.}
}
\value{
A \code{\link{GRanges-class}} object.
}
\description{
Load structural variant calls
}
\author{
David Porubsky
}
