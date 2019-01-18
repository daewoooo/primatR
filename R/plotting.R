#' Plot categorized metacolumns of ranged data.
#' 
#' This function counts categorcal variables stored as metacolumn in \code{\link{GRanges-class}} object.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param colName A metacolumn name to be sumarized and plotted.
#' @param facedID A metacolumn name to be used to split data in sub-plots.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

plotColumnCounts <- function(gr, colName='gen', facetID=NULL) {
  plt.df <- as.data.frame(gr)
  ## Use user defined name of the column to plot
  if (is.null(facetID)) {
    data.tab <- plt.df %>% group_by(.dots=eval(colName)) %>% summarize(counts = length(eval(parse(text = colName))))
  } else {
    data.tab <- plt.df %>% group_by(.dots=c(eval(facetID), eval(colName) )) %>% summarize(counts = length(eval(parse(text = colName))))
  }
  
  plt <- ggplot(data.tab) + geom_col(aes_string(x=eval(colName), y='counts', fill=eval(colName)))
  plt <- plt + geom_text(data=data.tab, aes_string(x=eval(colName), y='counts', label='counts'), vjust=0)
  plt <- plt + scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"), guide='none')
  plt <- plt + theme_bw() + xlab("")
  
  if (!is.null(facetID)) {
    plt <- plt + facet_grid(eval(parse(text = facetID)) ~ .)
  }  
  
  return(plt)
}


#' Plot categorized metacolumns of ranged data per chromosome.
#' 
#' This function counts categorcal variables stored as metacolumn in \code{\link{GRanges-class}} object per each chromosome.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param colName A metacolumn name to be sumarized and plotted.
#' @param facedID A metacolumn name to be used to split data in sub-plots.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

plotColumnCountsPerChr <- function(gr, colName='gen', facetID=NULL, normChrSize=FALSE) {
  plt.df <- as.data.frame(gr)
  ## Use user defined name of the column to plot
  if (is.null(facetID)) {
    data.tab <- plt.df %>% group_by(.dots=c('seqnames', eval(colName) )) %>% summarize(counts = length(eval(parse(text = colName))))
    #data.tab <- plt.df %>% group_by(.dots='seqnames') %>% summarize(counts = n()) # counts per chr only
  } else {
    data.tab <- plt.df %>% group_by(.dots=c('seqnames', eval(facetID), eval(colName) )) %>% summarize(counts = length(eval(parse(text = colName))))
  }
  
  if (normChrSize) {
      seq.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
      data.tab <- data.tab %>% group_by(seqnames) %>% mutate(ChrSizeNorm=counts / seq.len[seqnames])
      
      plt <- ggplot(data.tab) + geom_col(aes_string(x='seqnames', y='ChrSizeNorm', fill=eval(colName)))
      plt <- plt + scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"))
      plt <- plt + theme_bw() + xlab("")
  } else {
      plt <- ggplot(data.tab) + geom_col(aes_string(x='seqnames', y='counts', fill=eval(colName)))
      plt <- plt + scale_fill_manual(values = brewer.pal(n = 9, name = "Set1"))
      plt <- plt + theme_bw() + xlab("")
  }
      
  if (!is.null(facetID)) {
    plt <- plt + facet_grid(eval(parse(text = facetID)) ~ .)
  }  
  
  return(plt)
}


#' Plot distribution of total sizes of all inversions per genotype. 
#'
#' This function plots sum of widths of all ranged data (total bp) per genotype and per chromosome.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

basesPerGenotypePerChr <- function(gr, normChrSize=FALSE) {
  plt.df <- as.data.frame(gr)
  data.tab <- plt.df %>% group_by(.dots=c("ID","seqnames","gen")) %>% summarize(counts = length(gen), inv.bases=sum(width))
  
  if (normChrSize) {
      seq.len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
      data.tab <- data.tab %>% group_by(seqnames) %>% mutate(ChrSizeNorm=inv.bases / seq.len[seqnames])
      plt <- ggplot(data.tab) + geom_col(aes(x=seqnames, y=ChrSizeNorm, fill=gen))
  } else {
      plt <- ggplot(data.tab) + geom_col(aes(x=seqnames, y=inv.bases, fill=gen))
  }    
      
  plt <- plt + scale_fill_manual(values = brewer.pal(n = 4, name = "Set1")[3:4], name="")
  plt <- plt + facet_grid(ID ~ .) + theme_bw()
  plt <- plt + xlab("Chromosomes") + ylab("Inverted bases (bp)")
  return(plt)
} 

#' Plot sorted size distribution of ranged data.
#' 
#' This function take \code{\link{GRanges-class}} object and plot width distribution of all ranges per variable stored in metacolumn field 'ID'
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param plotUniqueBases Set to \code{TRUE} if you want to plot size of unique bases per range.
#' @param violin Set to \code{TRUE} if you want plot violo_plot instead of dot_plot.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

rangesSizeDistribution <- function(gr, plotUniqueBases=FALSE, violin=FALSE) {
  inv.sizes.ord <- order(width(gr), decreasing = FALSE)
  if (plotUniqueBases) {
    size.dist.df <- data.frame(x=1:length(inv.sizes.ord), size=width(gr)[inv.sizes.ord], uniqueBases=gr$TotalUniqueBases[inv.sizes.ord], ID=gr$ID[inv.sizes.ord])
  } else {
    size.dist.df <- data.frame(x=1:length(inv.sizes.ord), size=width(gr)[inv.sizes.ord], ID=gr$ID[inv.sizes.ord])
  }  
  
  if (violin) {
    plt <- ggplot(size.dist.df) + geom_violin(aes(x=ID, y=size, fill=ID), trim = FALSE)
    plt <- plt + geom_dotplot(aes(x=ID, y=size), binaxis='y', stackdir='center', dotsize=0.05, binwidth = 1)
    plt <- plt + scale_fill_manual(values = brewer.pal(n = 4, name = "Set1"), name="")
  } else {
    plt <- ggplot(size.dist.df) + geom_point(aes(x=x, y=size, color=ID))
    if (plotUniqueBases) {
      plt <- plt + geom_point(aes(x=x, y=uniqueBases), color='gray')
    }
  }
  
  plt <- plt + scale_y_continuous(breaks=c(1000,10000,100000,1000000), labels = comma, trans = 'log10')
  plt <- plt + geom_hline(yintercept = c(10000, 1000000), linetype="dashed")
  plt <- plt + scale_color_manual(values = brewer.pal(n = 4, name = "Set1"), name="")
  plt <- plt + xlab("Size sorted inversions") + ylab("Inversion size (log10)")
  return(plt)
}  


#' Plots scatter of event counts to chromosome size
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

eventsPerChrSizeScatter <- function(gr, bsgenome) {
  if (any(is.na(seqlengths(gr))) & is.null(bsgenome)) {
    message("Chromosome lengths are missing. Please submit BSgenome object of the genome you want to plot.")
  }
  
  ## Load BSgenome
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      suppressPackageStartupMessages(library(bsgenome, character.only=T))
      bsgenome <- as.object(bsgenome) # replacing string by object
    }
  }
  
  plt.df <- as.data.frame(gr)
  seq.len <- seqlengths(bsgenome)[seqlevels(gr)]
  data.tab <- plt.df %>% group_by(.dots='seqnames') %>% summarize(counts = n()) %>% mutate(ChrLen=seq.len[seqnames])

  data.tab %>% ggplot() + geom_point(aes(x=ChrLen, y=counts)) + 
  geom_text(data=data.tab, aes(x=ChrLen, y=counts, label=seqnames), vjust=-0.5, hjust=-0.1) +
  scale_x_continuous(labels = comma) +  
  xlab("Chromosome size (bp)")
}


#' Plot genome-wide distribution of ranged data.
#' Ranges are color by the 'ID' column.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

genomewideRangesIdeo <- function(gr, userTrack=NULL, userTrackGeom='rect', bsgenome=NULL) {
  if (any(is.na(seqlengths(gr))) & is.null(bsgenome)) {
    message("Chromosome lengths are missing. Please submit BSgenome object of the genome you want to plot.")
  }
  
  ## Load BSgenome
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      suppressPackageStartupMessages(library(bsgenome, character.only=T))
      bsgenome <- as.object(bsgenome) # replacing string by object
    }
  }
  
  seq.len <- seqlengths(bsgenome)[paste0('chr', c(1:22, 'X'))]
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=paste0('chr', c(1:22, 'X')))
  
  gr$level <- GenomicRanges::disjointBins(gr)
  plt.df <- as.data.frame(gr)
  
  plt <- ggplot(ideo.df) + geom_linerange(aes(x=-1, ymin=0, ymax=length), size=2, color="black")
  if (!is.null(userTrack)) {
    strand(userTrack) <- "*"
    userTrack <- reduce(userTrack)
    ## Remove seqlevels not present in the data or different than chr1:22 and X
    chroms2keep <- seqlevels(userTrack)[seqlevels(userTrack) %in% paste0('chr', c(1:22, 'X'))]
    userTrack <- keepSeqlevels(userTrack, chroms2keep, pruning.mode = 'coarse')
    userTrack.df <- as(userTrack, 'data.frame')
    if (userTrackGeom == 'rect') {
      plt <- plt + geom_linerange(data=userTrack.df, aes(x=-1, ymin=as.numeric(start), ymax=as.numeric(end)), size=2, color="darkgoldenrod1")
    } else if (userTrackGeom == 'point') {
      userTrack.df$midpoint <- userTrack.df$start + ((userTrack.df$end - userTrack.df$start)/2)
      plt <- plt + geom_point(data=userTrack.df, aes(x=-1, y=midpoint), size=1, color="darkgoldenrod1")
    }
  }
  plt <- plt + coord_flip() + facet_grid(seqnames ~ ., switch = 'y')
  plt <- plt + geom_linerange(data=plt.df , aes(x=level, ymin=start, ymax=end+250000, color=ID), size=1)
  plt <- plt + scale_color_manual(values =  RColorBrewer::brewer.pal(n = 4, name = "Set1"), name="")
  plt <- plt + scale_y_continuous(expand = c(0,0))
  plt <- plt + theme_void()
  plt <- plt + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  plt <- plt + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plt <- plt + theme(strip.text.y = element_text(angle = 180))
  return(plt)
}


#' Plot size distribution of ranges in all Venn partitions.
#'
#' @param venn A \code{data.frame} object as an output of \pkg{VennDiagram} package.
#' @param overlaps A \code{\link{GRanges-class}} object with overlaps determined by the same idx column ID.
#' @return A \code{ggplot} object.
#' @author David Porubsky

plotVennPartitions <- function(venn=NULL, overlaps=NULL) {
  partitions <- list()
  for (i in 1:nrow(venn)) {
    partition <- venn[i,]
    partition.idx <- unlist(partition$..values..)
    partition.ranges <- overlaps[overlaps$sub.group %in% partition.idx]
    partition.ranges$partitionID <- partition$..set..
    partitions[[i]] <- as.data.frame(partition.ranges)
  }
  
  plt.df <- do.call(rbind, partitions)
  plt <- ggplot(plt.df) + geom_boxplot(aes(x=ID, y=width, fill=ID))
  plt <- plt + scale_fill_manual(values = brewer.pal(n = 3, name = "Set1"), name="")
  plt <- plt + scale_y_continuous(breaks=c(1000,10000,100000,1000000), labels = comma, trans = 'log10')
  plt <- plt + geom_hline(yintercept = c(10000, 100000, 1000000), linetype="dashed")
  plt <- plt + xlab("") + ylab("Inversion size (log10)")
  plt <- plt + facet_grid(partitionID ~ .)
  return(plt)
}


#' Plot genotype distribution of ranges in all Venn partitions.
#'
#' @param venn A \code{data.frame} object as an output of \pkg{VennDiagram} package.
#' @param overlaps A \code{\link{GRanges-class}} object with overlaps determined by the same idx column ID.
#' @return A \code{ggplot} object.
#' @author David Porubsky

plotVennPartitionsGen <- function(venn=NULL, overlaps=NULL) {
  partitions <- list()
  for (i in 1:nrow(venn)) {
    partition <- venn[i,]
    partition.idx <- unlist(partition$..values..)
    partition.ranges <- overlaps[overlaps$idx %in% partition.idx]
    partition.ranges$partitionID <- partition$..set..
    partitions[[i]] <- as.data.frame(partition.ranges)
  }
  plt.df <- do.call(rbind, partitions)
  
  data.tab <- plt.df %>% group_by(.dots=c('partitionID', 'ID', 'gen')) %>% summarize(counts = length(gen))
  plt <- ggplot(data.tab, aes_string(x='ID', y='counts', group='gen')) 
  plt <- plt + geom_col(aes_string(fill='gen'), position = 'dodge')
  plt <- plt + geom_text(aes_string(label='counts'), vjust=0,  position=position_dodge(width = 1))
  plt <- plt + scale_fill_manual(values = brewer.pal(n = 3, name = "Set1"), name="")
  plt <- plt + xlab("")
  plt <- plt + facet_grid(partitionID ~ .)
  return(plt)
}


#' Plot number of overlaps between query and subject ranges.
#'
#' @param query.gr A \code{\link{GRanges-class}} object.
#' @param subject.gr A \code{\link{GRanges-class}} object. 
#' @param facedID A metacolumn name to be used to split data in sub-plots.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

plotOverlapWithRanges <- function(query.gr=NULL, subject.gr=NULL, facetID=NULL) {
  counts <- countOverlaps(query.gr, subject.gr)
  query.gr$userTrackOverlapCount <- counts
  
  inv.sizes.ord <- order(width(query.gr), decreasing=TRUE)
  
  plt.df <- as.data.frame(query.gr[inv.sizes.ord])
  plt.df$x <- 1:nrow(plt.df) 
  
  plt <- ggplot(plt.df) + geom_col(aes_string(x='x', y='width'), color='black')
  plt <- plt + geom_col(aes_string(x='x', y='userTrackOverlapCount'), color='red', inherit.aes=FALSE)
  plt <- plt + scale_y_continuous(breaks=c(1000,10000,100000,1000000), labels = comma, trans = 'log10')
  plt <- plt + geom_hline(yintercept = c(10000, 1000000), linetype="dashed")
  plt <- plt + xlab("Size sorted inversions (bp)") + ylab("Gene count | Inversion size")
  
  return(plt)
} 


#' Plot genome-wide distribution of plus and minus reads.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param bin.len A size of bins to count directional reads in.
#' @param colors A colors for plus and minus reads.
#' @return A \code{ggplot} object.
#' @importFrom gtools mixedsort
#' @importFrom scales comma
#' @author David Porubsky
#' @export

plotCompositeIdeo <- function(gr, bin.len = 200000, colors = c('#EFEDF5', '#68228B')) {
  ## Sort reads by position and by chromosome
  gr <- sort(gr, ignore.strand=TRUE)
  seqlevels(gr) <- gtools::mixedsort(seqlevels(gr))
  #seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[seqlevels(gr)]
  
  chrom.lengths <- seqlengths(gr)
  
  ## Bin the data
  chrom.lengths.floor <- floor(chrom.lengths/ bin.len) * bin.len
  binned.data <- unlist(tileGenome(chrom.lengths.floor, tilewidth = bin.len), use.names = FALSE)
  seqlengths(binned.data) <- chrom.lengths
  
  Watsonreads <- GenomicRanges::countOverlaps(binned.data, gr[strand(gr)=='-']) 
  Crickreads <- GenomicRanges::countOverlaps(binned.data, gr[strand(gr)=='+'])
  bothreads <- Watsonreads + Crickreads
  
  mcols(binned.data)$bothreads <- bothreads
  mcols(binned.data)$Watsonreads <- Watsonreads
  mcols(binned.data)$Crickreads <- Crickreads
  
  dfplot.reads <- as.data.frame(binned.data)
  
  Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.99)
  Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.99)
  ## Set outlier bins to the limit
  dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
  dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
  
  my.theme <- theme(legend.position="none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    strip.text.y = element_text(angle = 180) ) 
  
  dfplot.reads$midpoint <- dfplot.reads$start + ( (dfplot.reads$end - dfplot.reads$start) %/% 2 )
  
  ## Construct ggplot
  dfplot.reads$mCrickreads <- -dfplot.reads$Crickreads
  ggplt <- ggplot(dfplot.reads)
  ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color=colors[1], size=1)
  ggplt <- ggplt + facet_grid(seqnames ~ ., switch = "y", scales = 'free')
  ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color=colors[2], size=1)
  ggplt <- ggplt +
    xlab("Genomic position") +
    ylab("Read counts") +
    scale_x_continuous(expand = c(0,0), labels = scales::comma) + my.theme
  return(ggplt)
}


#' Plot 'dotplot' of two sequences.
#'
#' This function specifically plots matched region of two sequences based on nucmer results.
#' 
#' @param nucmer.coords A coordinates from nucmer output. [RUN: nucmer --coords ...] 
#' @param genome.coords Set to \code{TRUE} if you want to work in genomic coordinates.
#' @param highlight.pos A set of postions to be highlighted on the x-axis.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' @export
plotNucmerCoords <- function(nucmer.coords = NULL, genome.coord = TRUE, highlight.pos = NULL, title = NULL, sd.track = NULL) {
  
  ## Helper function
  remapCoord <- function(x = NULL, new.range = NULL) {
    offset <- (min(new.range) + x[1]) - 1
    dist <- cumsum(diff(x))
    new.x <- c(offset, (offset + dist))
    return(new.x)
  } 
  
  ## Read in coordinates from nucmer output
  coords <- read.table(nucmer.coords, skip=5, stringsAsFactors = FALSE)
  coords.df <- data.frame(s1.start=coords$V1,
                          s1.end=coords$V2,
                          s2.start=coords$V4,
                          s2.end=coords$V5,
                          s1.id=coords$V12,
                          s2.id=coords$V13, 
                          stringsAsFactors = FALSE)
  
  ## Translate sequence coordinates to genome-wide coordinates
  if (genome.coord) {
    s1.region <- unique(coords.df$s1.id)
    s2.region <- unique(coords.df$s2.id)
    s1.chr <- strsplit(s1.region, ":")[[1]][1]
    s1.region <- strsplit(s1.region, ":")[[1]][2]
    s2.region <- strsplit(s2.region, ":")[[1]][2]
    s1.range <- as.numeric( strsplit(s1.region, "-")[[1]] )
    s2.range <- as.numeric( strsplit(s2.region, "-")[[1]] )
    coords.df$s1.start <- remapCoord(x = coords.df$s1.start, new.range = s1.range)
    coords.df$s1.end <- remapCoord(x = coords.df$s1.end, new.range = s1.range)
    coords.df$s2.start <- remapCoord(x = coords.df$s2.start, new.range = s2.range)
    coords.df$s2.end <- remapCoord(x = coords.df$s2.end, new.range = s2.range)
  }
  
  ## Color segments by directionality
  coords.df$dir <- 'rev'
  forw.mask <- (coords.df$s1.start < coords.df$s1.end) & (coords.df$s2.start < coords.df$s2.end)
  coords.df$dir[forw.mask] <- 'forw'
  
  ## Plot segments
  plt <- ggplot2::ggplot(coords.df, aes(x=s1.start,xend=s1.end,y=s2.start,yend=s2.end, color=dir)) + 
    geom_segment() +
    xlab(unique(coords.df$s1.id)) +
    ylab(unique(coords.df$s2.id)) +
    theme_bw() + 
    theme(aspect.ratio=1) + #force x and y axis to have the same proportion
    scale_color_manual(values = c('chartreuse4', 'darkgoldenrod2'))
  ## Highlight user defined postion on x axis
  if (!is.null(highlight.pos)) {
    plt <- plt + geom_vline(xintercept = highlight.pos, color='black', linetype = 2)
  }
  ## Highlight regions of SDs on x axis
  if (!is.null(sd.track) & genome.coord) {
    s1.gr <- GRanges(seqnames=s1.chr, ranges=IRanges(start=s1.range[1], end=s1.range[2]))
    sd.track <- subsetByOverlaps(sd.track, s1.gr)
    sd.track.df <- as.data.frame(sd.track)
    # make sure that sd regions do not extend over x axis region
    sd.track.df$start[sd.track.df$start < start(s1.gr)] <- start(s1.gr)
    sd.track.df$end[sd.track.df$end > end(s1.gr)] <- end(s1.gr)
    
    plt <- plt + geom_rect(data=sd.track.df  ,aes(xmin=start, xmax=end, ymin=Inf, ymax=-Inf), fill='gray', alpha=0.25, inherit.aes = FALSE)
  }
  ## Highlight user defined postion on x axis
  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  
  return(plt)
}


#' Plot HIC contact matrix.
#'
#' This function imports HIC reads and plot contact matrix based on user defined resolution. 
#' 
#' @param resolution A coordinates from nucmer output. [RUN: nucmer --coords ...] 
#' @param bsgenome A set of postions to be highlighted on the x-axis.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @inheritParams bam2GRanges
#' @inheritParams plotNucmerCoords
#' @author David Porubsky
#' @export
plotHICregional <- function(bamfile, region = NULL, min.mapq = 10, resolution = 50000, highlight.pos = NULL, bsgenome = NULL) {
  
  ## Read in paired-end reads
  message("Loading reads for region: ", as.character(region), " ...")
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, 
                                                                       index=paste0(bamfile,'.bai'), 
                                                                       param=Rsamtools::ScanBamParam(tag="XA", 
                                                                                                     which=range(region), 
                                                                                                     what=c('mapq','flag'),
                                                                                                     flag=Rsamtools::scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE)
                                                                       )
  )
  )
  ## Convert data into GRanges objects
  data.first <- as(GenomicAlignments::first(data.raw), 'GRanges') #mate1
  data.last <- as(GenomicAlignments::last(data.raw), 'GRanges') #mate2
  
  ## Filter reads that spans user-defined region
  hits.first <- findOverlaps(query = data.first, subject = region)
  hits.last <- findOverlaps(query = data.last, subject = region)
  mask <- intersect(queryHits(hits.first), queryHits(hits.last))
  data.first <- data.first[mask]
  data.last <- data.last[mask]
  
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    mask <- data.first$mapq >= min.mapq & data.last$mapq >= min.mapq
    data.first <- data.first[mask]
    data.last <- data.last[mask]
  }
  
  plots <- list()
  for (res in resolution) {
    message("    Constructing contact matrix and plotting for binsize ", res, "bp ...")
    ## Make genomic bins
    bins <- makeBins(bsgenome = bsgenome, chromosomes = as.character(seqnames(region)), binsize = res, stepsize = res)
    bins <- subsetByOverlaps(bins, region)
    
    ## Preparing all possible combination of genomic bins
    pairs <- t(combn(length(bins), 2))
    
    ## Get overlaps between reads and genomic bins
    hits.first <- findOverlaps(subject = bins, query = data.first) #mate1
    hits.last <- findOverlaps(subject = data.last, query = bins) #mate2
    
    ## Split reads indices by genomic bins indices
    frags.idx.per.bin.first <- split(queryHits(hits.first), subjectHits(hits.first)) #mate1
    frags.idx.per.bin.last <- split(subjectHits(hits.last), queryHits(hits.last)) #mate2
    
    ## Count shared mate1 and mate2 indices for all combinations of genomic bins
    link.counts1 <- mapply(M1=frags.idx.per.bin.first[pairs[,1]], M2=frags.idx.per.bin.last[pairs[,2]], function(M1, M2) length(intersect(M1, M2)))
    link.counts2 <- mapply(M1=frags.idx.per.bin.last[pairs[,1]], M2=frags.idx.per.bin.first[pairs[,2]], function(M1, M2) length(intersect(M1, M2)))
    link.counts <- link.counts1 + link.counts2
    
    ## Construct data.frame with counts per contact bins
    link.df <- data.frame('M1'=as.data.frame(bins[pairs[,1]])[,1:3], 'M2'=as.data.frame(bins[pairs[,2]])[,1:3], value=link.counts)
    
    ## Plot binned contacts
    plt <- ggplot(link.df[link.df$value > 0,]) + 
      geom_rect(aes(xmin=M1.start, xmax=M1.end, ymin=M2.start, ymax=M2.end, fill=value))
    if (!is.null(highlight.pos)) {
      plt <- plt + geom_vline(xintercept = highlight.pos, color='red') +
        geom_hline(yintercept = highlight.pos, color='red')
    }
    plt <- plt + scale_fill_gradient(low = "white", high = "black", trans = 'log', name="Links (log)") +
      scale_x_continuous(labels = comma) +
      scale_y_continuous(labels = comma) +
      xlab("Genomic position (bp)") +
      ylab("Genomic position (bp)") +
      annotate(geom="text", x=max(link.df$M1.start), y=min(link.df$M1.start), label=paste0("Resolution: ", res, "bp"), hjust=1) +
      theme(aspect.ratio=1) #Make sure plot is always a square
    plots[[length(plots) + 1]] <- plt 
  }
  message("DONE!!!")
  return(plots)
}
