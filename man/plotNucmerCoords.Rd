% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotNucmerCoords}
\alias{plotNucmerCoords}
\title{Plot 'dotplot' of two sequences.}
\usage{
plotNucmerCoords(
  nucmer.coords = NULL,
  genome.coord = TRUE,
  highlight.pos = NULL,
  title = NULL,
  sd.track = NULL,
  shape = "segm"
)
}
\arguments{
\item{nucmer.coords}{A coordinates from nucmer output. [RUN: nucmer --coords ...]}

\item{highlight.pos}{A set of postions to be highlighted on the x-axis.}

\item{title}{A character string to use as a title for the plot.}

\item{sd.track}{A segmental dulication track to be highlight at the plot.}

\item{shape}{A shape used to plot aligned sequences: Either 'segm' or 'point'.}

\item{genome.coords}{Set to \code{TRUE} if you want to work in genomic coordinates.}
}
\value{
A \code{\link[ggplot2:ggplot]{ggplot}} object.
}
\description{
This function specifically plots matched region of two sequences based on nucmer results.
}
\author{
David Porubsky
}
