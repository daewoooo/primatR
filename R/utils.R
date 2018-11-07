#' Insert chromosome for in case it's missing
#' 
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} object with an additional metadata column containing chromosome name with 'chr'. 
insertchr <- function(gr) {
  mask <- which(!grepl('chr', seqnames(gr)))
  mcols(gr)$chromosome <- as.character(seqnames(gr))
  mcols(gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(gr)$chromosome[mask])
  mcols(gr)$chromosome <- as.factor(mcols(gr)$chromosome)
  return(gr)
}


as.object <- function(x) {
  return(eval(parse(text=x)))
}

#' Resize genomic ranges to user defined multiple of their original size
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param times A number of how many times to extend each range by its own width.
#' @param bsgenome A reference genome to get lengths of standard chromosomes (1-22 adn X).
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @export
resizeRanges <- function(gr, times=2, bsgenome) {
  ## Load BSgenome
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      suppressPackageStartupMessages(library(bsgenome, character.only=T))
      bsgenome <- as.object(bsgenome) # replacing string by object
    }
  }
  
  seqlengths(gr) <- seqlengths(bsgenome)[seqlevels(gr)]
  
  for (i in seq_along(gr)) {
    range <- gr[i]
    range.seqLen <- seqlengths(range)[unique(seqnames(range))]
    extension <- width(range)*times
    new.start <- max(start(range) - extension, 1)
    new.end <- min(end(range) + extension, range.seqLen)
    start(gr[i]) <- new.start
    end(gr[i]) <- new.end
  }  
  return(gr)
}


#' This function calculates coverage per set of sequencing fragments or reports disjoint position per fragment.
#'
#' @param grl A \code{\link{GRangesList-class}} object.
#' @param ID A unique identifier to add as a metacolumn in returned \code{\link{GRangesList}} object.
#' @param coverage A logical value (TRUE|FALSE) if to report collapsed coverage or single reads.
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @export
coveragePerRegion <- function(grl, ID="", coverage=TRUE) {
  region.covs <- GRangesList()
  for (i in seq_along(grl)) {
    gr <- grl[[i]]
    if (coverage) {
      cov.plus <- Biostrings::coverage(gr[strand(gr) == '+'])
      cov.minus <- Biostrings::coverage(gr[strand(gr) == '-'])
    
      cov.plus.ranges <- as(cov.plus[unique(seqnames(gr))], 'GRanges')
      cov.minus.ranges <- as(cov.minus[unique(seqnames(gr))], 'GRanges')
      cov.plus.ranges <- cov.plus.ranges[-c(1, length(cov.plus.ranges))]
      cov.minus.ranges <- cov.minus.ranges[-c(1, length(cov.minus.ranges))]
      cov.minus.ranges$score <- cov.minus.ranges$score * -1
      strand(cov.plus.ranges) <- "+"
      strand(cov.minus.ranges) <- "-"
    
      cov.ranges <- c(cov.plus.ranges, cov.minus.ranges)
      cov.ranges$ID <- unique(gr$ID)
      cov.ranges$region.ID <- unique(gr$region.ID)
    } else {
      plus.ranges <- gr[strand(gr) == "+"]
      minus.ranges <- gr[strand(gr) == "-"]
      plus.ranges$level <- disjointBins(plus.ranges)
      minus.ranges$level <- -disjointBins(minus.ranges)
      cov.ranges <- c(plus.ranges, minus.ranges)
    }
    region.covs[[i]] <-  cov.ranges
  }
  return(region.covs)
}


# Reduce overlapping ranges
collapseOverlaps <- function(gr) {
  reduced.gr <- reduce(gr)
  if (ncol(mcols(gr))>0) {
    mcols(reduced.gr) <- mcols(gr)[length(reduced.gr),]
  }
  return(reduced.gr)
}