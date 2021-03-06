% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeRanges.R
\name{removeRanges}
\alias{removeRanges}
\title{Function to subtract set of genomic ranges from other genomic ranges}
\usage{
removeRanges(gr = NULL, remove.gr = NULL, collapse.gaps = 0, min.region = 0)
}
\arguments{
\item{gr}{A \code{\link{GRanges-class}} object of ranges from which remove.gr are removed.}

\item{remove.gr}{A \code{\link{GRanges-class}} object of ranges to subtract from gr.}

\item{collapse.gaps}{Size of the gaps between ranges to be merged.}

\item{min.region}{Minimum size regions to be reported.}
}
\value{
A \code{\link{GRanges-class}} object.
}
\description{
Function to subtract set of genomic ranges from other genomic ranges
}
\author{
David Porubsky
}
