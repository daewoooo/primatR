#' Export aligned reads into a bedgraph
#'
#' This function takes as an input aligned reads in BAM and exports them into a UCSC formated bedgraph.
#'
#' @param bamfile Bamfile with aligned reads.
#' @param outputdirectory A directory to save constructed bedgraphs.
#' @param mapq Minimum mapping quality when importing from BAM files.
#' @param filt.flag Filter out reads with given flag.
#' @param min.read.len Minimum length of mapped reads to import.
#' @param blacklist A \code{\link{GRanges-class}} object of regions to exclude reads from.
#' @return \code{NULL}
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author David Porubsky
#' @export 

exportBedGraph <- function(bamfile, outputdirectory=".", mapq=10, filt.flag=0, min.read.len=5000, blacklist=NULL) {
  filename <- basename(bamfile)
  filename <- gsub(filename, pattern = "\\.bam$", replacement = "")
  
  ## Load reads from BAM
  data.raw <- GenomicAlignments::readGAlignments(bamfile, param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'qname'), flag=scanBamFlag(isDuplicate=FALSE)))
  data <- as(data.raw, 'GRanges')

  ## Filter by mapq
  if (mapq > 0) {
    data <- data[data$mapq >= mapq]
  }  
  ## Filter by flag
  if (filt.flag > 0) {
    bit.flag <- bitwAnd(filt.flag, data$flag)
    mask <- bit.flag == 0
    data <- data[mask]
  }  
  ## filter by read length
  data <- data[width(data) >= min.read.len]
  ## Get breakpoint ID
  data$ID <- sapply(data$qname, function(x) strsplit(x, "__")[[1]][2])
  ## Split by breakpoint
  data.grl <- GenomicRanges::split(data, data$ID)
  
  ## Set RGB colors
  #colors <- RColorBrewer::brewer.pal(n = length(data.grl), name = 'Set1')
  #color.rgb <- grDevices::col2rgb(colors)
  #color.rgb <- apply(color.rgb, 2, function(x) paste(x, collapse = ","))
  
  ## Set random RGB colors
  color.rgb <- sapply(1:length(data.grl), function(x) round(runif(3, 0, 255)))
  color.rgb <- apply(color.rgb, 2, function(x) paste(x, collapse = ","))
  
  ## Create file to save results
  savefile <- paste0(filename, "_", ".bedgraph.gz")
  destination <- file.path(outputdirectory, savefile)
  savefile.gz <- gzfile(destination, 'w')
  for (i in seq_along(data.grl)) {
    gr <- data.grl[[i]]
    gr <- GenomeInfoDb::keepSeqlevels(gr, paste0('chr', c(1:22,'X')), pruning.mode = 'coarse')
    cov <- GenomicRanges::coverage(gr)
    cov.gr <- as(cov, 'GRanges')
    
    ## Filter out regions that overlap with user defined blacklisted regions
    if (class(blacklist) == "GRanges") {
      hits <- GenomicRanges::findOverlaps(cov.gr, blacklist)
      cov.gr <- cov.gr[-queryHits(hits)]
    }
    
    ## Get 2 max covered region for each breakpoint
    cov.gr.perChr <- split(cov.gr, seqnames(cov.gr))
    bases.total <- sapply(cov.gr.perChr, function(x) sum(x$score))
    bases.total <- sort(bases.total, decreasing = T)[1:2]
    max.bases <- paste0(names(bases.total), "-", bases.total)
    max.bases <- paste(max.bases, collapse = "_")
    
    header<- paste0('track type=bedGraph name=', unique(gr$ID), ' description=SplittedPacBioReads_', filename, '_', max.bases,' visibility=full color=', color.rgb[i])
    write.table(header, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
    
    bedF <- data.frame(cov.gr)[c('seqnames','start','end','score')]
    write.table(bedF, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
  }
  close(savefile.gz)
}

exportBedGraphPairs <- function(bamfile, outputdirectory=".", mapq=10, filt.flag=0, min.read.len=5000, blacklist=NULL) {
  filename <- basename(bamfile)
  filename <- gsub(filename, pattern = "\\.bam$", replacement = "")
  
  ## Load paired reads from BAM
  data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'qname'), flag=scanBamFlag(isDuplicate=FALSE))) 
  data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
  data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
  
  ## Filter by mapq
  if (mapq > 0) {
    mask <- data.first$mapq >= mapq &  data.last$mapq >= mapq
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
  ## Filter by read length
  mask <- width(data.first) >= min.read.len & width(data.last) >= min.read.len
  data.first <- data.first[mask]
  data.last <- data.last[mask]
  ## Filter out blacklisted regions
  if (!is.null(blacklist)) {
    hits.first <- findOverlaps(query = data.first, subject = blacklist)
    hits.last <- findOverlaps(query = data.last, subject = blacklist)
    filt <- c(queryHits(hits.first), queryHits(hits.last))
    if (length(filt) > 0) {
      data.first <- data.first[-filt]
      data.last <- data.last[-filt]
    }
  }  
  
  ## Create final data object
  data <- data.first
  data$mate <- data.last
  
  ## Get breakpoint ID
  data$ID <- sapply(data$qname, function(x) strsplit(x, "__")[[1]][2])
  ## Split by breakpoint
  data.grl <- GenomicRanges::split(data, data$ID)
  
  ## Set random RGB colors
  color.rgb <- sapply(1:length(data.grl), function(x) round(runif(3, 0, 255)))
  color.rgb <- apply(color.rgb, 2, function(x) paste(x, collapse = ","))
  
  ## Create file to save results
  savefile <- paste0(filename, "_", ".bedgraph.gz")
  destination <- file.path(outputdirectory, savefile)
  savefile.gz <- gzfile(destination, 'w')
  for (i in seq_along(data.grl)) {
    gr <- data.grl[[i]]
    gr.mate1 <- gr[,0]
    gr.mate2 <- gr$mate[,0]
    ## Keep only standard chromosomes
    gr.mate1 <- GenomeInfoDb::keepSeqlevels(gr.mate1, paste0('chr', c(1:22,'X')), pruning.mode = 'coarse')
    gr.mate2 <- GenomeInfoDb::keepSeqlevels(gr.mate2, paste0('chr', c(1:22,'X')), pruning.mode = 'coarse')
    ## Calculate coverage separately for mate1 and mate2
    cov.mate1 <- GenomicRanges::coverage(gr.mate1)
    cov.mate1.gr <- as(cov.mate1, 'GRanges')
    cov.mate2 <- GenomicRanges::coverage(gr.mate2)
    cov.mate2.gr <- as(cov.mate2, 'GRanges')
    
    ## Get 2 max covered region for each breakpoint
    cov.gr <- c(cov.mate1.gr, cov.mate2.gr)
    cov.gr.perChr <- split(cov.gr, seqnames(cov.gr))
    bases.total <- sapply(cov.gr.perChr, function(x) sum(x$score))
    bases.total <- sort(bases.total, decreasing = T)[1:2]
    max.bases <- paste0(names(bases.total), "-", bases.total)
    max.bases <- paste(max.bases, collapse = "_")
    
    ## Print mate1 to bedgraph
    header<- paste0('track type=bedGraph name=', unique(gr$ID), '_mate1 description=SplittedPacBioReads_', filename, '_', max.bases,' visibility=full color=', color.rgb[i])
    write.table(header, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
    
    bedF <- data.frame(cov.mate1.gr)[c('seqnames','start','end','score')]
    write.table(bedF, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
    
    ## Print mate2 to bedgraph
    header<- paste0('track type=bedGraph name=', unique(gr$ID), '_mate2 description=SplittedPacBioReads_', filename, '_', max.bases,' visibility=full color=', color.rgb[i])
    write.table(header, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
    
    bedF <- data.frame(cov.mate2.gr)[c('seqnames','start','end','score')]
    write.table(bedF, file=savefile.gz, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
  }
  close(savefile.gz)
}
