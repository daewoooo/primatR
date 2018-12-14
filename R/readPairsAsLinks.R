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

processLinks <- function(inputfolder = ".", min.mapq = 60, filt.flag = 3328, min.reads = 10, chromosomes = NULL, bsgenome = NULL, blacklist = NULL) {

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
    
    ## Remove positions overlapping with SDs
    hits.first <- IRanges::findOverlaps(data.first, seg.dup.gr)
    hits.last <- IRanges::findOverlaps(data.last, seg.dup.gr)
    mask <- unique(GenomicRanges::sort(c(queryHits(hits.first), queryHits(hits.last))))
    data.first <- data.first[-mask]
    data.last <- data.last[-mask]
    
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
    data <- data.first
    data$to.gr <- data.last
    data$ID <- sapply(data$qname, function(x) strsplit(x, "__")[[1]][2])
    fragments[[i]] <- data
  }
  
  ## Merge all data for plotting
  data.plt <- unlist(fragments, use.names = FALSE)
  assessed.breaks <- length(unique(data.plt$ID))
  
  ## TODO filter random links here???
  binned.gr <- makeBins(bsgenome = bsgenome, chromosomes = chromosomes, binsize = 100000, stepsize = 50000)
  from.gr <- data.plt[,'qname']
  to.gr <- data.plt$to.gr[,'qname']
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
  data.plt <- data.plt[filt.idx]
  return(list(links = data.plt, assessed.breaks = assessed.breaks))
}  


#' Filter and plot significant connections between read pairs
#' 
#' This function takes condidate links and filter all spurious links. Also this function exports table of links 
#' as well as plot links onto the genome reference.
#'
#' @param gr.links Minimum number of reads to support a link.
#' @inheritParams processLinks
#' @author David Porubsky
#' @export

plotLinks <- function(gr.links, min.reads = 10, chromosomes = NULL) {

  ## Get number of assessed breaks from the link object
  assessed.breaks <- gr.links$assessed.breaks
  
  ## Prepare data.frame for plotting
  plt.df <- as.data.frame(gr.links$links)
  plt.df$break.chr <- sapply(plt.df$ID, function(x) gsub(strsplit(x, "-")[[1]][1], pattern = 'chr', replacement = ''))
  
  ## Add x_coord by selecting larger number being the start of the link
  plt.df$x <- pmax(plt.df$start, plt.df$to.gr.start)
  plt.df$xend <- pmin(plt.df$start, plt.df$to.gr.start)
  ## Remove links with the same start and end position
  plt.df$link <- paste0(plt.df$x, "_", plt.df$xend)
  plt.df <- plt.df[!duplicated(plt.df$link),]
  ## Remove links with the same start and end position
  plt.df <- plt.df[plt.df$x != plt.df$xend,]
  
  ## Filter only user defined chromosomes
  mask <- plt.df$seqnames %in% chromosomes & plt.df$to.gr.seqnames %in% chromosomes
  plt.df <- plt.df[mask,]
  
  ## Impose chromosome levels to observed data
  plt.df$y <- as.character(plt.df$seqnames)
  plt.df$yend <- as.character(plt.df$to.gr.seqnames)
  plt.df$y <- match(plt.df$y, as.character(chromosomes))
  plt.df$yend <- match(plt.df$yend, as.character(chromosomes))
  
  ## Make sure that at least one break resides on chromosome of origin
  #mask <- plt.df$break.chr != plt.df$y & plt.df$break.chr != plt.df$yend
  plt.df <- plt.df %>% filter(break.chr == y | break.chr == yend)
  
  ## Add links between chromosomes
  chr.from <- pmax(plt.df$y, plt.df$yend)
  chr.to <- pmin(plt.df$y, plt.df$yend)
  plt.df$chr.link <- paste0(chr.from, "_", chr.to)
  
  ## Get number of unique reads per link
  qnames.list <- split(plt.df$qname, plt.df$chr.link)
  qnames.counts.perLink <- sapply(qnames.list , countUniqueReadIDs)
  ## Export links
  export.link.table <- plt.df
  export.link.table <- export.link.table[!duplicated(export.link.table$chr.link),]
  export.link.table$link.unique.reads <-  qnames.counts.perLink
  ## Filter links with minimal unique read support
  signif.links <- names(qnames.counts.perLink)[qnames.counts.perLink >= min.reads]
  plt.df.filt <- plt.df[plt.df$chr.link %in% signif.links,]
  
  ## Get number of intra- versus inter-chromosomal links
  total.breaks <- length(unique(plt.df.filt$ID))
  breaks.df <- split(plt.df.filt, plt.df.filt$ID)
  intra.links <- 0
  inter.links <- 0
  for (i in seq_along(breaks.df)) {
    break.df <- breaks.df[[i]]
    links <- unique(break.df$chr.link)
    for (j in seq_along(links)) {
      link <- links[j]
      chroms <- strsplit(link, "_")[[1]]
      
      if (chroms[1] == chroms[2]) {
        intra.links <- intra.links + 1
      } else {
        inter.links <- inter.links + 1
      }
    }  
  }
  
  ## Prepare genome-wide ideogram
  seq.len <- seqlengths(bsgenome)[chromosomes]
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chromosomes)
  ideo.df$levels <- 1:length(seq.len)
  ideo <- ggplot(ideo.df) + geom_linerange(aes(x=levels, ymin=0, ymax=length), size=1, color='black')
  
  ## Add links to the ideogram
  plt <- ideo +
    geom_curve(data=plt.df.filt, aes(x = y, y = start, xend = yend, yend = to.gr.start), color="red", curvature = 0.25) + 
    scale_y_continuous(labels = comma) +
    scale_x_continuous(breaks = ideo.df$levels, labels = ideo.df$seqnames) +
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    theme_bw() +
    ylab("Genomic position (bp)") +
    xlab("Chromosome")
  
  ## Add some statistics
  plt <- plt + annotate(geom="text", x=23, y=250000000, label=paste0("Assessed breaks: ", assessed.breaks, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-10000000, label=paste0("Intra-chromosomal breaks: ", intra.links, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-20000000, label=paste0("Inter-chromosomal breaks: ", inter.links, "\n"), color="black", hjust=1)
  
  return(list(plot=plt, links=export.link.table))
}  