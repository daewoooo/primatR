#' Calculate genotype in a set of genomic ranges
#'
#' This function exports genotype information for a set of user defined genomic ranges.
#'
#' @param regions A \code{\link{GRanges-class}} object that contains genomic regions to be genotyped.
#' @param directional.reads A \code{\link{GRanges-class}} object that contains directional Strand-seq reads.
#' @param blacklist A \code{\link{GRanges-class}} object containing regions to be excluded from genotyping.
#' @param index User defined character string to be appended to calculated genotypes.
#' @inheritParams getGenotype
#' @return A \code{\link{GRanges-class}} object with appended count of plus and minus reads and the most likely genoptype.
#' @author David Porubsky
#' @export
#'
genotypeRegions <- function(regions=NULL, directional.reads=NULL, blacklist=NULL, min.reads=5, alpha=0.05, index=NULL) {
  ## Initialize output GRanges object
  regions.genot <- regions
  mcols.df <- data.frame(plus.reads = rep(0, length(regions.genot)),
                         minus.reads = rep(0, length(regions.genot)),
                         genoT = rep('lowReads', length(regions.genot)),
                         stringsAsFactors = FALSE)
  ## Subset directional reads by overlaps with region of interest
  frags <- IRanges::subsetByOverlaps(directional.reads, regions)
  ## Remove fragments overlapping with blackilisted regions
  if (!is.null(blacklist)) {
    blacklist <- GenomicRanges::reduce(blacklist)
    frags <- IRanges::subsetByOverlaps(frags, blacklist, invert = TRUE)
  }
  ## Split fragments by the region they overlap with
  #hits <- IRanges::findOverlaps(frags, regions, select = 'first')
  #frags.grl <- GenomicRanges::split(frags, hits)
  hits <- IRanges::findOverlaps(frags, regions)
  frags.grl <- GenomicRanges::split(frags[queryHits(hits)], subjectHits(hits))
  ## Get genotype for each region
  genoT.l <- lapply(frags.grl, function(x) getGenotype(gr = x, min.reads = min.reads, alpha = alpha))
  genoT.df <- do.call(rbind, genoT.l)
  ## Add regions that could have been genotype to output data.frame
  mcols.df[rownames(genoT.df),] <- genoT.df
  ## Append user defined index to column names
  if (!is.null(index)) {
      if (nchar(index) > 0) {
        colnames(mcols.df) <- paste0(c('plus.reads', 'minus.reads', 'genoT'), "_", index)
      }  
  }  
  ## Genotypes to the output GRanges object
  mcols(regions.genot) <- c(mcols(regions.genot), mcols.df)
  return(regions.genot)
}

#' Calculate genotype for a specific genomic region
#'
#' Based on directional Strand-seq reads from a genomic region reports most likely genotype based on binomial probabilites of allowed strand states (WW, CC or WC).
#'
#' @param gr A \code{\link{GRanges-class}} object that contains directional Strand-seq reads.
#' @param min.reads Minimal number of reads to pursue genotyping.
#' @param alpha Estimated level of background in Strand-seq reads.
#' @return A \code{data.frame} of binomial probabilities for a given counts of plus and minus reads for a single cell. (rows=reads/genomic segments, cols=strand states)
#' @author David Porubsky
#' @export
#' 
getGenotype <- function(gr, min.reads=5, alpha=0.05) {
  if (length(gr) >= min.reads) {
    ## Counts plus and minus reads
    plus.reads <- length(gr[strand(gr) == '+'])
    minus.reads <- length(gr[strand(gr) == '-'])
    ## Calculate binomial probabilities
    bin.probs <- primatR::countProb(minusCounts = minus.reads, plusCounts = plus.reads, alpha = alpha, log = TRUE)
    max.prob <- which.max(bin.probs)
    if (max.prob == 1) {
      gen <- 'HOM'
    } else if (max.prob == 2) {
      gen <- 'REF'
    } else {
      gen <- 'HET'
    }
    ## Calculate ww ratio to determine genotype
    # ww.ratio <- (plus.reads - minus.reads) / (plus.reads + minus.reads)
    # if (ww.ratio < -0.5) {
    #   gen <- 'HOM'
    # } else if (ww.ratio > 0.5) {
    #   gen <- 'REF'
    # } else {
    #   gen <- 'HET'
    # }
  } else {
    ## Counts plus and minus reads
    plus.reads <- length(gr[strand(gr) == '+'])
    minus.reads <- length(gr[strand(gr) == '-'])
    gen <- "lowReads"
  } 
  return(data.frame(plus.reads = plus.reads, minus.reads = minus.reads, genoT = gen, stringsAsFactors = FALSE))
}


#' Calculate probability of a genomic region having WW, CC or WC state.
#'
#' Exports most likely strand state (WW, CC or WC) based on counts of plus and minus reads and allowed level of background 'alpha'.
#'
#' @param minusCounts Minus (Watson) read counts aligned to PacBio reads.
#' @param plusCounts Plus (Crick) read counts aligned to PacBio reads.
#' @param alpha Estimated level of background in Strand-seq reads.
#' @return A \code{matrix} of binomial probabilities for a given counts of plus and minus reads for a single cell. (rows=reads/genomic segments, cols=strand states)
#' @author David Porubsky, Maryam Ghareghani
#' @export
#' 
countProb <- function(minusCounts, plusCounts, alpha=0.1, log=FALSE) {
  
  sumCounts <- minusCounts + plusCounts
  #calculate that given region is WW
  prob.ww <- stats::dbinom(minusCounts, size = sumCounts, prob = 1-alpha, log = log)
  
  #calculate that given region is CC
  prob.cc <- stats::dbinom(minusCounts, size = sumCounts, prob = alpha, log = log)
  
  #calculate that given region is WC
  prob.mix <- stats::dbinom(minusCounts, size = sumCounts, prob = 0.5, log = log)
  
  prob.m <- cbind(prob.ww, prob.cc, prob.mix)
  
  return(prob.m)
}
