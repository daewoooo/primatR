% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateCN.R
\name{estimateCN}
\alias{estimateCN}
\title{Estimate copy number (CN) in genomic regions}
\usage{
estimateCN(frags, regions, zlim = 3.291)
}
\arguments{
\item{frags}{A \code{\link{GRanges-class}} object that contains directional Strand-seq reads.}

\item{regions}{A \code{\link{GRanges-class}} object containing regions defined by breakpoint boundaries.}

\item{zlim}{The number of stdev that the deltaW must pass the peakTh (ensures only significantly higher peaks are considered).}
}
\value{
A \code{\link{GRanges-class}} object containing submitted regions with an extra metacolumn \code{CN.estim}.
}
\description{
This function calculates ratio between observed and expected number of reads in localized genomic regions as CN estimates.
}
\author{
David Porubsky
}
