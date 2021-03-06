% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getReciprocalOverlaps.R
\name{getReciprocalOverlaps}
\alias{getReciprocalOverlaps}
\title{Function to calculate percentage of reciprocal overlap between two sets of ranges.}
\usage{
getReciprocalOverlaps(query, subject, thresh = 50, report = "both", index = "")
}
\arguments{
\item{query}{A \code{\link{GRanges-class}} object.}

\item{subject}{A \code{\link{GRanges-class}} object.}

\item{thresh}{A percentage threshold for a required overlap.}

\item{report}{Set to 'both' if you want to return merged query and subject ranges or select only one of them.}

\item{index}{An user defiend ID to disntiguish calculated overlaps.}
}
\value{
A \code{\link{GRanges-class}} object.
}
\description{
Function to calculate percentage of reciprocal overlap between two sets of ranges.
}
\author{
David Porubsky
}
