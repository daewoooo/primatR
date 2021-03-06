% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{eventsPerChrSizeScatter}
\alias{eventsPerChrSizeScatter}
\title{Plots scatter of event counts to chromosome size}
\usage{
eventsPerChrSizeScatter(gr, bsgenome, colBy = NULL, facetID = NULL, lm = FALSE)
}
\arguments{
\item{gr}{A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.}

\item{bsgenome}{A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.}

\item{colBy}{A metacolumn name to be used to color data by.}

\item{facetID}{A metacolumn name to be used to split data in sub-plots.}

\item{lm}{Set to TRUE if regression line should be added to the plot.}
}
\value{
A \code{ggplot} object.
}
\description{
Plots scatter of event counts to chromosome size
}
\author{
David Porubsky
}
