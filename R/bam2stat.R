#' Get basic BAM file statistics 
#'
#' This function takes as an input aligned reads in BAM and exports some useful statistics like depth of coverage
#' and percentage of genome being coverged by reads.
#'
#' @param bamfile Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param chunkSize Set for big BAMs to process them in smaller chunks.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param filt.flag Filter out reads with a given flag.
#' @param filt.alt Set to \code{TRUE} if you want to filter out alternative alignments defined in 'XA' tag.
#' @return \code{data.frame}
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @author David Porubsky
#' @export 

bam2stat <- function(bamfile, bamindex=paste0(bamfile, '.bai'), chromosomes=NULL, chunkSize=NULL, min.mapq=10, filt.flag=3328, filt.alt=FALSE) {
  filename <- basename(bamfile)
  message("Reading BAM: ", filename)
  ## Check if bamindex exists
  bamindex.raw <- sub('\\.bai$', '', bamindex)
  bamindex <- paste0(bamindex.raw,'.bai')
  if (!file.exists(bamindex)) {
    bamindex.own <- Rsamtools::indexBam(bamfile)
    warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
    bamindex <- bamindex.own
  }
  file.header <- Rsamtools::scanBamHeader(bamfile)[[1]]
  chrom.lengths <- file.header$targets
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  if (length(chroms2use) == 0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff) > 0) {
    diffs <- paste0(diff, collapse=', ')
    warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
  }
  ## Chromosome to load
  gr <- GenomicRanges::GRanges(seqnames=chroms2use, ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
  
  if (!is.null(chunkSize)) {
    message("Processing BAM in chunks ...")
    bam.in <- Rsamtools::BamFile(bamfile, yieldSize=chunkSize)
    open(bam.in)
    
    n <- 0 # number of reads examined
    ## read in and count chunks of data, until done
    total.bases <- 0
    covered.pos <- 0
    total.reads <- 0
    repeat {
      ## Load the reads in chunks
      data.raw <- GenomicAlignments::readGAlignments(bam.in, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what=c('mapq', 'flag')))# input next chunk
      
      if (length(data.raw) == 0) {break}
      
      n <- n + length(data.raw)
      message("    Processed reads: ", n)
      
      data <- as(data.raw, "GRanges")
      ## Filter by mapping quality
      if (min.mapq > 0) {
        data <- data[data$mapq >= min.mapq]
      }
      ## Filter by flag
      if (filt.flag > 0) {
        bit.flag <- bitwAnd(filt.flag, data$flag)
        mask <- bit.flag == 0
        data <- data[mask]
      } 
      ## Filter out reads with alternative alignments
      if (filt.alt) {
        data <- data[is.na(data$XA)]
      }
      
      cov.rle <- GenomicRanges::coverage(data)
      cov.gr <- as(cov.rle, "GRanges")
      
      ## Summary stat ##
      total.bases <- total.bases + sum(cov.gr$score)
      covered.pos <- covered.pos + sum( width(cov.gr[cov.gr$score > 0]) )
      total.reads <- total.reads + length(data)
    }
    ## Create final data.frame
    summary.df <- data.frame(filename=filename, total.bases=total.bases,covered.pos=covered.pos, total.reads=total.reads)
    ## finish and return result
    close(bam.in)
  
  } else {
    ## Read in all reads
    data.raw <- GenomicAlignments::readGAlignments(bamfile, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what=c('mapq', 'flag')))
    data <- as(data.raw, "GRanges")
    
    ## Filter by mapping quality
    if (min.mapq > 0) {
      data <- data[data$mapq >= min.mapq]
    }
    ## Filter by flag
    if (filt.flag > 0) {
      bit.flag <- bitwAnd(filt.flag, data$flag)
      mask <- bit.flag == 0
      data <- data[mask]
    } 
    ## Filter out reads with alternative alignments
    if (filt.alt) {
      data <- data[is.na(data$XA)]
    }
    
    cov.rle <- GenomicRanges::coverage(data)
    cov.gr <- as(cov.rle, "GRanges")
    
    ## Summary stat ##
    covered.ranges <- cov.gr[cov.gr$score > 0]
    total.aligned.bases <- sum(width(covered.ranges) * covered.ranges$score)
    covered.genomic.pos <- sum(width(covered.ranges))
    total.reads <- length(data)
    median.read.len <- median(width(data))
    summary.df <- data.frame(filename=filename, 
                             total.aligned.bases=total.aligned.bases,
                             covered.genomic.pos=covered.genomic.pos, 
                             total.reads=total.reads,
                             median.read.len=median.read.len)
  }
  return(summary.df)
}