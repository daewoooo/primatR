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

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} with two additional metadata columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
  cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
  cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
  names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
  gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
  return(gr)
}

## Helper function
as.object <- function(x) {
  return(eval(parse(text=x)))
}

#' Resize genomic ranges to user defined multiple of their original size
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param times A number of how many times to extend each range by its own width.
#' @param bsgenome A reference genome to get lengths of standard chromosomes (1-22 adn X).
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @author David Porubsky
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
#' @author David Porubsky
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
    
      cov.ranges <- sort(c(cov.plus.ranges, cov.minus.ranges), ignore.strand=TRUE)
      cov.ranges$ID <- unique(gr$ID)
      cov.ranges$region.ID <- unique(gr$region.ID)
    } else {
      plus.ranges <- gr[strand(gr) == "+"]
      minus.ranges <- gr[strand(gr) == "-"]
      plus.ranges$level <- disjointBins(plus.ranges)
      minus.ranges$level <- -disjointBins(minus.ranges)
      cov.ranges <- sort(c(plus.ranges, minus.ranges), ignore.strand=TRUE)
    }
    region.covs[[i]] <-  cov.ranges
  }
  return(region.covs)
}


#' Split genome into bins
#' 
#' This function splid genome into a equaly-sized user defined genomic intervals.
#'
#' @param bsgenome A reference genome to get lengths of genomic sequences (eg. GRCh38).
#' @param fai A FASTA index to get lengths of genomic sequences.
#' @param chromosomes A user defined set of chromosomes for binning (eg. 'chr1')
#' @param binsize A size of the genomic bin to split genome into.
#' @param stepsize A size of the genomic interval to move each bin. For non-overlapping bins use the same size as binsize.
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges. 
#' @author David Porubsky
#' @export
makeBins <- function(bsgenome=NULL, fai=NULL, chromosomes, binsize=100000, stepsize=binsize/2) {
  
  if (!is.null(bsgenome)) {
    chr.lengths <- seqlengths(bsgenome)[chromosomes]
  } else if (!is.null(fai)) {
    ## Get contigs/scaffolds names and sizes from fasta index
    fai.tab <- utils::read.table(fai)
    fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
    chr.lengths <- fai.tab$V2
    names(chr.lengths) <- fai.tab$V1
    chr.lengths <- chr.lengths[names(chr.lengths) %in% chromosomes]
  } else {
    warning("Please submit chromosome lengths in a form of BSgenome object or fasta index (.fai)!!!")
  } 
  
  bins <- GRangesList()
  seqlevels(bins) <- chromosomes
  for (i in seq_along(chr.lengths)) {
    chr.len <- chr.lengths[i]
    
    bin.starts <- seq(from = 1, to = chr.len-binsize, by = stepsize)
    bin.ends <- seq(from = binsize, to = chr.len, by = stepsize)
    
    chr.bins <- GRanges(seqnames=names(chr.len), ranges=IRanges(start=bin.starts, end=bin.ends))
    bins[[i]] <- chr.bins
  }
  bins <- unlist(bins, use.names = FALSE)
  seqlengths(bins) <- chr.lengths[chromosomes]
  return(bins)
}


#' Reduce set of overlapping genomic ranges
#'
#' @param gr A \code{\link{GRanges-class}} object of set of genomic intervals.
#' @return A \code{\link{GRanges-class}} object with resized original set of ranges.
#' @author David Porubsky
#' @export
collapseOverlaps <- function(gr) {
  reduced.gr <- GenomicRanges::reduce(gr)
  if (ncol(mcols(gr)) > 0) {
    mcols(reduced.gr) <- unique(mcols(gr)) #[length(reduced.gr),]
  }
  return(reduced.gr)
}


#' Count number of unique read IDs in a set of reads.
#'
#' @param readIDs Set of read IDs with appened extra info using '__' as a delimiter.
#' @author David Porubsky
#' @export
countUniqueReadIDs <- function(readIDs) {
  read.id <- sapply(readIDs, function(x) strsplit(x, "__")[[1]][1])
  unique.ids <- unique(read.id)
  return(length(unique.ids))
} 


#' Collapses consecutive set of ranges with the same value
#' 
#' @param gr A \code{\link{GRanges}} object.
#' @param id.field A field column used to collapse ranges with the same value.
#' @param measure.field A field column that contains measured value to be sum up.
#' @author David Porubsky
#' @export
#' 
collapseBins <- function(gr, id.field=0, measure.field=NULL) {
  ## Include seqnames into the unique ID
  unique.ID <- paste0(seqnames(gr), '_', mcols(gr)[,id.field])
  ## Get continous runs of the same unique ID
  unique.ID.runLength <- runLength(Rle(unique.ID))
  ind.last <- cumsum(unique.ID.runLength) ##get indices of last range in a consecutive(RLE) run of the same value
  ind.first <- c(1,cumsum(unique.ID.runLength) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
  collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
  names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
  ## Sum values in user defined measure fields(s)
  if (!is.null(measure.field)) {
    run.ID  <- factor(rep(seq_along(unique.ID.runLength), unique.ID.runLength))
    for (field in measure.field) {
      mcols(collapsed.gr)[,field] <-  as.numeric(tapply(mcols(gr)[,field], run.ID, sum))
      #collapsed.gr$C <-  tapply(gr$C, run.ID, sum)
      #collapsed.gr$W <-  tapply(gr$W, run.ID, sum)
    }
  }
  return(collapsed.gr)
}


#' Bind consecutive pairs values in vector into a matrix with two columns
#' 
#' @param v A \code{\link{vector}} object.
#' 
#' @author David Porubsky
#' @export
#'
reformat <- function(v) {
  lines <- list()
  for (i in 1:(length(v)-1)) {
    lines[[length(lines)+1]] <- c(v[i], v[i+1])
  }
  return(do.call(rbind, lines))
}


#' Expand genomic Ranges into a set of region boundaries
#' 
#' @param gr A \code{\link{GRanges}} object.
#' @param sort If set to \code{TRUE} reported \code{\link{GRanges}} object will be sorted by position.
#' @author David Porubsky
#' @export
#'
getRegionBoundaries <- function(gr, sort=TRUE) {
  ## Extract start and end position of each genomic range
  starts.gr <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start=start(gr), end=start(gr)))
  ends.gr <- GRanges(seqnames = seqnames(gr), ranges = IRanges(start=end(gr), end=end(gr)))
  ## Region ID to both start and end position
  starts.gr$region.ID <- as.character(gr)
  ends.gr$region.ID <- as.character(gr)
  ## Sort and export region boundaries
  if (sort) {
    all.pos.gr <- sort(c(starts.gr, ends.gr), ignore.strand=TRUE)
  } else {
    starts.gr$ord <- 1:length(starts.gr)
    ends.gr$ord <- 1:length(ends.gr)
    all.pos.gr <- c(starts.gr, ends.gr)
    all.pos.gr <- all.pos.gr[order(all.pos.gr$ord), 1]
  }  
  return(all.pos.gr)
}

#' Make sure that directionality of set of ranges reflect desired majority strand
#'
#' @param ranges A \code{\link{GRanges}} object.
#' @param majority.strand A desired majority strand directionality to be reported.
#' @param strand.only If set to \code{TRUE} simple character string of strand directionality ('+', '-') is reported.
#' @author David Porubsky
#' @export
#'
syncRangesDir <- function(ranges, majority.strand = '+', strand.only = FALSE) {
  ## Define majority and minority strand
  if (majority.strand == '+') {
    minority.strand = '-'
  } else if (majority.strand == '-') {
    minority.strand = '+'
  } else {
    stop("Parameter 'majority.strand' can only be defined as '+' or '-' !!!")
  }
  ## Check if the input object is GRangesList or single GRanges object
  if (grepl(class(ranges), pattern = 'GRangesList')) {
    for (i in seq_along(ranges)) {
      gr <- ranges[[i]]
      ## Flip directionality based to make sure majority strand covers the most bases
      if (sum(width(gr[strand(gr) == majority.strand])) > sum(width(gr[strand(gr) == minority.strand]))) {
        ranges[[i]] <- gr
      } else {
        gr.new <- gr
        strand(gr.new)[strand(gr) == majority.strand] <- minority.strand
        strand(gr.new)[strand(gr) == minority.strand] <- majority.strand
        ranges[[i]] <- gr.new
      }
    }
    if (strand.only) {
      return(as.character(strand(unlist(ranges, use.names = FALSE))))
    } else {
      return(ranges)
    }
  } else if (class(ranges) == 'GRanges') {
    gr <- ranges
    ## Flip directionality based to make sure majority strand covers the most bases
    if (sum(width(gr[strand(gr) == majority.strand])) > sum(width(gr[strand(gr) == minority.strand]))) {
      if (strand.only) {
        return(as.character(strand(gr)))
      } else {
        return(gr)
      }  
    } else {
      gr.new <- gr
      strand(gr.new)[strand(gr) == majority.strand] <- minority.strand
      strand(gr.new)[strand(gr) == minority.strand] <- majority.strand
      if (strand.only) {
        return(as.character(strand(gr.new)))
      } else {
        return(gr.new)
      }
    }
  } else {
    stop("Only valid 'GRanges' or 'GRangesList' object can be processed !!!")
  }
}
