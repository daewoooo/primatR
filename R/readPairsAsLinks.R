#' Get significant connections between read pairs
#' 
#' This function takes split read mappings of PacBio reads and keeps only genomic regions with at least 10 unique 
#' reads in a given (100kb) genomic bin.
#'
#' @param inputfolder A directory that contains set of BAM files to be processed.
#' @param min.reads Minimum number of unique reads to support a link.
#' @param blacklist An \code{\link{GRanges-class}} object that contains regions to be filtered out.
#' @return A \code{ggplot} object.
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @inheritParams bam2stat
#' @inheritParams makeBins
#' @author David Porubsky
#' @export

readPairsAsLinks <- function(inputfolder = ".", min.mapq = 60, filt.flag = 3328, min.reads = 10, chromosomes = NULL, bsgenome = NULL, blacklist = NULL) {

  ## Get list of BAM files in the inputfolder
  bamfiles <- list.files(inputfolder, pattern = "\\.bam$", full.names = TRUE)
  
  fragments <- GenomicRanges::GRangesList()
  for (i in seq_along(bamfiles)) {
    bam <- bamfiles[[i]]
    filename <- basename(bam)
    message("Processing BAM: ", filename)
    
    ## Load read pairs
    suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bam, param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'qname'))) )
    ## Convert to GRanges
    data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
    data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
    
    ## Remove reads overlapping with blacklisted regions
    hits.first <- IRanges::findOverlaps(data.first, seg.dup.gr)
    hits.last <- IRanges::findOverlaps(data.last, seg.dup.gr)
    mask <- unique(GenomicRanges::sort(c(queryHits(hits.first), queryHits(hits.last))))
    if (length(mask) > 0) {
      data.first <- data.first[-mask]
      data.last <- data.last[-mask]
    }  
    
    ## Filter by mapq
    if (min.mapq > 0) {
      mask <- data.first$mapq >= min.mapq & data.last$mapq >= min.mapq
      data.first <- data.first[mask]
      data.last <- data.last[mask]
    }
    
    ## Filter by flag
    if (filt.flag > 0) {
      bit.flag.first <- bitwAnd(filt.flag, data.first$flag)
      bit.flag.last <- bitwAnd(filt.flag, data.last$flag)
      mask <- bit.flag.first == 0 & bit.flag.last == 0
      data.first <- data.first[mask]
      data.last <- data.last[mask]
    } 
    
    ## Export fragments
    if (length(data.first) > 0) {
      data <- data.first
      data$to.gr <- data.last
      data$ID <- sapply(data$qname, function(x) strsplit(x, "__")[[1]][2])
      data$bam.name <- filename
      fragments[[length(fragments) + 1]] <- data
    }  
  }
  
  ## Merge all data
  fragments.all <- unlist(fragments, use.names = FALSE)
  assessed.breaks <- length(unique(fragments.all$ID))
  assessed.bams <- unique(fragments.all$bam.name)
  
  ## Count number of reads in 100kb long bins in order to filter out noisy regions
  binned.gr <- makeBins(bsgenome = bsgenome, chromosomes = chromosomes, binsize = 100000, stepsize = 50000)
  from.gr <- fragments.all[,'qname']
  to.gr <- fragments.all$to.gr[,'qname']
  all.gr <- c(from.gr, to.gr)
  
  hits <- IRanges::findOverlaps(binned.gr, all.gr)
  readIDs.per.overlap <- split(all.gr$qname[subjectHits(hits)], queryHits(hits))
  unique.readIDs.per.overlap <- sapply(readIDs.per.overlap, countUniqueReadIDs)
  ## Get bins with minimal number of unique reads
  bin.idx <- as.numeric( names(unique.readIDs.per.overlap)[unique.readIDs.per.overlap >= min.reads] )
  select.regions <- binned.gr[bin.idx]
  ## Collapse overlapping ranges
  filt.regions <- GenomicRanges::reduce(select.regions)
  
  ## Find overlapping ranges
  from.gr.hits <- IRanges::findOverlaps(from.gr, filt.regions)
  to.gr.hits <- IRanges::findOverlaps(to.gr, filt.regions)
  filt.idx <- IRanges::intersect(queryHits(from.gr.hits), queryHits(to.gr.hits))
  fragments.all <- fragments.all[filt.idx]
  return(list(links = fragments.all, assessed.breaks = assessed.breaks, assessed.bams = assessed.bams))
}