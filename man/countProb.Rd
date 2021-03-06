% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genotypeRegions.R
\name{countProb}
\alias{countProb}
\title{Calculate probability of a genomic region having WW, CC or WC state.}
\usage{
countProb(minusCounts, plusCounts, alpha = 0.1, log = FALSE)
}
\arguments{
\item{minusCounts}{Minus (Watson) read counts aligned to PacBio reads.}

\item{plusCounts}{Plus (Crick) read counts aligned to PacBio reads.}

\item{alpha}{Estimated level of background in Strand-seq reads.}
}
\value{
A \code{matrix} of binomial probabilities for a given counts of plus and minus reads for a single cell. (rows=reads/genomic segments, cols=strand states)
}
\description{
Exports most likely strand state (WW, CC or WC) based on counts of plus and minus reads and allowed level of background 'alpha'.
}
\author{
David Porubsky, Maryam Ghareghani
}
