#' Plot categorized metacolumns of ranged data.
#' 
#' This function counts categorcal variables stored as metacolumn in \code{\link{GRanges-class}} object.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param colName A metacolumn name to be sumarized and plotted.
#' @param facedID A metacolumn name to be used to split data in sub-plots.
#' @param colors A user defined set of colors used for plotting.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

plotColumnCounts <- function(gr, colName='gen', facetID=NULL, colors=NULL) {
  plt.df <- as.data.frame(gr)
  ## Use user defined name of the column to plot
  if (is.null(facetID)) {
    data.tab <- plt.df %>% group_by(.dots=eval(colName)) %>% summarize(counts = length(eval(parse(text = colName))))
  } else {
    data.tab <- plt.df %>% group_by(.dots=c(eval(facetID), eval(colName) )) %>% summarize(counts = length(eval(parse(text = colName))))
  }
  
  ## Set colors
  col.categs <- length(unique(plt.df$ID))
  if (is.null(colors)) {
    col <- brewer.pal(n = col.categs, name = "Set1")
  } else if (length(colors) < col.categs) {
    message("Not enough colors specified!!! Please provide ", col.categs, "distict colors")
  } else {
    col <- colors
  }
  
  plt <- ggplot(data.tab) + geom_col(aes_string(x=eval(colName), y='counts', fill=eval(facetID)))
  plt <- plt + geom_text(data=data.tab, aes_string(x=eval(colName), y='counts', label='counts'), vjust=0)
  plt <- plt + scale_fill_manual(values = col, guide='none')
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
#' @param colors A user defined set of colors used for plotting.
#' @param title A user defined title of the plot.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

rangesSizeDistribution <- function(gr, plotUniqueBases=FALSE, violin=FALSE, colors=NULL, title=NULL) {
  inv.sizes.ord <- order(width(gr), decreasing = FALSE)
  if (plotUniqueBases) {
    size.dist.df <- data.frame(x=1:length(inv.sizes.ord), size=width(gr)[inv.sizes.ord], uniqueBases=gr$TotalUniqueBases[inv.sizes.ord], ID=gr$ID[inv.sizes.ord], stringsAsFactors = FALSE)
  } else {
    size.dist.df <- data.frame(x=1:length(inv.sizes.ord), size=width(gr)[inv.sizes.ord], ID=gr$ID[inv.sizes.ord], stringsAsFactors = FALSE)
  }  
  
  ## Get median size distribution per ID
  size.dist.df <- size.dist.df %>% group_by(ID) %>% mutate(med_size = median(size))
  
  ## Set colors
  col.categs <- length(unique(size.dist.df$ID))
  if (is.null(colors)) {
    col <- brewer.pal(n = col.categs, name = "Set1")
  } else if (length(colors) < col.categs) {
    message("Not enough colors specified!!! Please provide ", col.categs, "distict colors")
  } else {
    col <- colors
  }
  
  if (violin) {
    plt <- ggplot(size.dist.df) + geom_violin(aes(x=ID, y=size, fill=ID), trim = FALSE) +
           geom_dotplot(aes(x=ID, y=size), binaxis='y', stackdir='center', dotsize=0.05, binwidth = 1) +
           scale_fill_manual(values = col, name="", guide='none') +
           geom_text(aes(x = ID, y = med_size, label=med_size), color="white", vjust=-0.5) +
           xlab("") 
  } else {
    plt <- ggplot(size.dist.df) + geom_point(aes(x=x, y=size, color=ID)) +
           facet_grid(ID ~ .) +
           geom_hline(aes(yintercept = med_size, group=ID), color="black") +
           geom_text(aes(x = 1, y = med_size, label=med_size), color="black", vjust=-0.5) +
           scale_color_manual(values = col, name="", guide='none') +
           xlab("Size sorted inversions") 
    if (plotUniqueBases) {
      plt <- plt + geom_point(aes(x=x, y=uniqueBases), color='gray')
    }
  }
  
  plt <- plt + 
         scale_y_continuous(breaks=c(1000,10000,100000,1000000), labels = comma, trans = 'log10') +
         geom_hline(yintercept = c(10000, 1000000), linetype="dashed") +
         ylab("Inversion size (log10)")
  
  if (!is.null(title) & is.character(title)) {
    plt <- plt + ggtitle(title)
  }
  
  return(plt)
}  


#' Plots scatter of event counts to chromosome size
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @param colBy A metacolumn name to be used to color data by.
#' @param facetID A metacolumn name to be used to split data in sub-plots.
#' @param lm Set to TRUE if regression line should be added to the plot.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

eventsPerChrSizeScatter <- function(gr, bsgenome, colBy=NULL, facetID=NULL, lm=FALSE) {
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
  
  if (lm) {
    colBy=NULL
    facetID=NULL
    message("Parameters 'colBy' and 'facetID' cannot be used when 'lm=TRUE'")
    l.mod = lm(counts ~ ChrLen, data=data.tab)
    suppressWarnings( conf.level <- predict(l.mod, interval="prediction", level = 0.95) )
    data.tab <- data.tab %>% mutate(resid=abs(resid(l.mod)), fitted=fitted(l.mod))
    data.tab <- cbind(data.tab, conf.level)
    
    r.sq <- paste0('R^2=', round(summary(l.mod)$r.squared, 3))
    r.sq.df <- as.data.frame(r.sq)
    
    plt <- data.tab %>% ggplot() + geom_line(aes_string(x='ChrLen', y='fitted')) +
      geom_line(aes_string(x='ChrLen', y='lwr'), color = "red", linetype = "dashed") +
      geom_line(aes_string(x='ChrLen', y='upr'), color = "red", linetype = "dashed") +
      geom_point(aes_string(x='ChrLen', y='counts', color='resid')) +
      geom_text(aes_string(x='ChrLen', y='counts', label='seqnames'), vjust=-0.5, hjust=-0.1) +
      geom_text(data = r.sq.df, aes(x=-Inf, y=Inf, label=r.sq), inherit.aes = F, hjust=-0.5, vjust=1) +
      scale_colour_gradient(low="blue", high="red") +
      scale_x_continuous(labels = comma) +
      labs(x="Chromosome size (bp)", y='counts', colour="Residuals")
  } else if (!is.null(colBy) & is.character(colBy)) {
    data.tab <- plt.df %>% group_by(.dots=c('seqnames', eval(colBy))) %>% summarize(counts = n()) %>% mutate(ChrLen=seq.len[seqnames])
    plt <- data.tab %>% ggplot(aes_string(x='ChrLen', y='counts', color=eval(colBy))) + 
      geom_point() + 
      geom_text(data=data.tab, aes_string(x='ChrLen', y='counts', label='seqnames'), vjust=-0.5, hjust=-0.1) +
      scale_x_continuous(labels = comma) +
      scale_color_manual(values = brewer.pal(n = 9, name = "Set1")) +
      xlab("Chromosome size (bp)")
  } else {
    plt <- data.tab %>% ggplot(aes_string(x='ChrLen', y='counts')) + 
      geom_point() + 
      geom_text(data=data.tab, aes_string(x='ChrLen', y='counts', label='seqnames'), vjust=-0.5, hjust=-0.1) +
      scale_x_continuous(labels = comma) +  
      xlab("Chromosome size (bp)")
  }
  ## Split the plot by facetID
  if (!is.null(facetID)) {
    plt <- plt + facet_grid(eval(parse(text = facetID)) ~ ., scales = 'free') #+ geom_smooth(method = "lm", se = FALSE)
  }  
  return(plt)
}


#' Plot genome-wide distribution of ranged data.
#' Ranges are color by the 'ID' column.
#'
#' @param gr A \code{\link{GRanges-class}} object with metadata columns to be summarized and plotted.
#' @param userTrack ...
#' @param userTrackGeom ...
#' @param colors A user defined set of colors used for plotting.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to provide chromosome lengths for plotting.
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export

genomewideRangesIdeo <- function(gr, userTrack=NULL, userTrackGeom='rect', colors=NULL, bsgenome=NULL) {
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
  
  ## Prepare ideo data for plotting
  seq.len <- seqlengths(bsgenome)[paste0('chr', c(1:22, 'X'))]
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=paste0('chr', c(1:22, 'X')))
  
  ## Prepare ranged data for plotting
  gr$level <- GenomicRanges::disjointBins(gr)
  plt.df <- as.data.frame(gr)
  
  ## Set colors
  col.categs <- length(unique(plt.df$ID))
  if (is.null(colors)) {
    col <- brewer.pal(n = col.categs, name = "Set1")
  } else if (length(colors) < col.categs) {
    message("Not enough colors specified!!! Please provide ", col.categs, "distict colors")
  } else {
    col <- colors
  }
  
  plt <- ggplot(ideo.df) + geom_linerange(aes(x=-1, ymin=0, ymax=length), size=2, color="gray36")
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
  plt <- plt + 
    coord_flip() + 
    facet_grid(seqnames ~ ., switch = 'y') +
    geom_linerange(data=plt.df , aes(x=level, ymin=start, ymax=end+250000, color=ID), size=1) +
    scale_color_manual(values =  col, name="") +
    scale_y_continuous(expand = c(0,0)) +
    theme_void() +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(strip.text.y = element_text(angle = 180))
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
#' @param title A character string to use as a title for the plot.
#' @param sd.track A segmental dulication track to be highlight at the plot.
#' @param shape A shape used to plot aligned sequences: Either 'segm' or 'point'.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' @export
plotNucmerCoords <- function(nucmer.coords = NULL, genome.coord = TRUE, highlight.pos = NULL, title = NULL, sd.track = NULL, shape='segm') {
  
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
  
  if (shape == 'segm') {
    ## Plot alignments
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.start,xend=s1.end,y=s2.start,yend=s2.end, color=dir)) + 
      geom_segment()
  } else if (shape == 'point') {
    coords.df$s1.midpoint <- coords.df$s1.start + ((coords.df$s1.end - coords.df$s1.start)/2)
    coords.df$s2.midpoint <- coords.df$s2.start + ((coords.df$s2.end - coords.df$s2.start)/2)
    plt <- ggplot2::ggplot(coords.df, aes(x=s1.midpoint, y=s2.midpoint, color=dir)) + 
      geom_point()
  }  
  plt <- plt +  
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
    sd.track <- reduce(sd.track)
    if (length(sd.track) > 0) {
      sd.track.df <- as.data.frame(sd.track)
      # make sure that sd regions do not extend over x axis region
      sd.track.df$start[sd.track.df$start < start(s1.gr)] <- start(s1.gr)
      sd.track.df$end[sd.track.df$end > end(s1.gr)] <- end(s1.gr)
      
      plt <- plt + geom_rect(data=sd.track.df  ,aes(xmin=start, xmax=end, ymin=Inf, ymax=-Inf), fill='gray', alpha=0.25, inherit.aes = FALSE)
    }  
  }
  ## Highlight user defined postion on x axis
  if (!is.null(title)) {
    plt <- plt + ggtitle(title)
  }
  return(plt)
}


#' Plot read-pair coverage profile per breakpoint
#'
#' This function takes as an input aligned read-pairs in BAM and plots coverage profile per breakpoint.
#' Breakpoint positions are stored as part of a readName.
#' e.g. readName__chr1-109148195-109148196__chr1_invDup_gorilla_roi4_1
#'
#' @param plot.ambig Bamfile with aligned reads.
#' @param plot.pairs Bamfile with aligned reads.
#' @param view.range Bamfile with aligned reads.
#' @inheritParams exportBedGraph 
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs first last getDumpedAlignments
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export 
plotReadPairCoverage <- function(bamfile, mapq=10, filt.flag=0, min.read.len=5000, plot.ambig=FALSE, plot.pairs=FALSE, blacklist=NULL, view.range=100000) {
  
  filename <- basename(bamfile)
  
  ## Load paired reads from BAM
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bamfile, param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'qname'), flag=scanBamFlag(isDuplicate=FALSE))) ) 
  data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
  data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
  
  ## Plot ambiguous read pairs if TRUE
  if (plot.ambig) {
    data.dumped <- GenomicAlignments::getDumpedAlignments()
    data.dumped.gr <- as(data.dumped, 'GRanges')
  }  
  
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
  breakpoints <- data.frame(break.id = unique(data$ID))
  ## Convert breakpoints to data.frame
  breakpoints.df <- tidyr::separate(breakpoints, break.id, sep = "-", into = c('chr','start','end'), convert = TRUE)
  
  ## Get region boundaries
  region.start <- min(breakpoints.df$start)
  region.end <- max(breakpoints.df$end)
  if (view.range > 0) {
    x.min <- region.start - view.range
    x.max <-region.end + view.range
  }
  
  ## Plot coverage profile per breakpoint
  data.grl <- split(data, data$ID)
  plt.data <- list()
  for (i in seq_along(data.grl)) {
    gr <- data.grl[[i]]
    ID <- unique(gr$ID)
    cov <- coverage(gr)
    cov.gr <- as(cov, 'GRanges')
    cov.gr$ID <- ID
    plt.data[[ID]] <- as.data.frame(cov.gr)
  }
  plt.df <- do.call(rbind, plt.data)
  
  ## Set random RGB colors
  #color.rgb <- sapply(1:length(data.grl), function(x) round(runif(3, 0, 255)))
  #color.rgb <- apply(color.rgb, 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  
  plots <- list()
  ## Prepare title
  title <- ggdraw() + draw_label(filename, fontface='bold')
  plots[[length(plots)+1]] <- title
  
  plt <- ggplot(plt.df) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=score, fill=ID)) + 
    geom_vline(xintercept = c(breakpoints.df$start)) +
    coord_cartesian(xlim = c(x.min, x.max)) +
    scale_fill_manual(values = brewer.pal(n = max(3, length(data.grl)), name = "Set1")) +
    theme(legend.position="bottom") +
    ggtitle("Breakpoint coverage")
  plots[[length(plots)+1]] <- plt
  
  if (plot.pairs) {
    plt.df <- as.data.frame(data)
    plt.df$level <- unlist(sapply(table(plt.df$ID), function(x) 1:x), use.names = FALSE)
    
    plt <- ggplot(plt.df) + 
      geom_segment(aes(x=start, xend=mate.start, y=level, yend=level), size=0.25) +
      geom_point(aes(x=start, y=level, color=strand), shape=108) +
      geom_point(aes(x=mate.start, y=level, color=mate.strand), shape=108) +
      geom_vline(xintercept = c(breakpoints.df$start)) +
      coord_cartesian(xlim = c(x.min, x.max)) + 
      theme(legend.position="bottom",
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle("Read-pair links") +
      xlab('')
    plots[[length(plots)+1]] <- plt
  }
  
  if (plot.ambig) {
    data.dumped.cov <- coverage(data.dumped.gr)
    data.dumped.cov.gr <- as(data.dumped.cov, "GRanges")
    plt.df <- as.data.frame(data.dumped.cov.gr)
    
    plt <- ggplot(plt.df) +
      geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=score, fill=ID), fill='orange') + 
      geom_vline(xintercept = c(breakpoints.df$start)) +
      coord_cartesian(xlim = c(x.min, x.max)) +
      theme(legend.position="bottom") +
      xlab('Genomic position (bp)') +
      ggtitle("Coverage of ambiguous read pairs")
    plots[[length(plots)+1]] <- plt
  }
  
  heights <- c(0.1, rep(1, length(plots)-1))
  
  final.plt <- cowplot::plot_grid(plotlist = plots, ncol = 1, rel_heights = heights)
  return(final.plt)
}


#' Plot BAM alignments around user defined genomic regions.
#'
#' @param bamfile Bamfile with aligned reads.
#' @param regions A \code{\link{GRanges-class}} object with genomic regions to process.
#' @param file A filename where final plot should be exported.
#' @return A \code{ggplot} object.
#' @importFrom dplyr group_by summarise "%>%"
#' @importFrom Rsamtools scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments cigarRangesAlongReferenceSpace cigarToRleList
#' @author David Porubsky
#' @export
#' 
plotAlignmentsPerRegion <- function(bamfile=NULL, regions=NULL, file=NULL) {
  ## Load BAM file
  message("Reading BAM file: ", basename(bamfile))
  data.raw <- GenomicAlignments::readGAlignments(bamfile, param=Rsamtools::ScanBamParam(what=c('cigar', 'mapq', 'flag', 'qname'), flag=scanBamFlag(isDuplicate=FALSE)))
  frags <- as(data.raw, 'GRanges')
  
  ## Go over user defined genomic location and plot BAM alignments at these regions
  plots <- list()
  for (i in seq_along(regions)) {
    invDup.region <- regions[i]
    regions.ID <- as.character(invDup.region)
    message("    Working on region: ", regions.ID)
    ## Select fragments from a genomic region
    frags.region <- IRanges::subsetByOverlaps(frags, invDup.region)
    if (length(frags.region) > 0) {
      ## Parse cigar string
      cigar.iranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(frags.region$cigar, flag = frags.region$flag, pos = start(frags.region))
      cigarRle <- GenomicAlignments::cigarToRleList(frags.region$cigar)
      ## Convert parsed cigar to GRanges
      cigar.gr <- GenomicRanges::GRanges(seqnames=as.character(seqnames(frags.region[1])), ranges=unlist(cigar.iranges))
      cigar.gr$cigar <- unlist(runValue(cigarRle))
      #qname.id <- sapply(frags.region$qname, function(x) strsplit(x, "\\.")[[1]][3])
      cigar.gr$qname <- factor(rep(frags.region$qname, lengths(cigar.iranges)), levels = unique(frags.region$qname))
      
      ## Restrict plotted region to the size of genomic regions of interest to right and left
      lookup.region <- primatR::resizeRanges(invDup.region, times = 1, bsgenome = BSgenome.Hsapiens.UCSC.hg38)
      cigar.gr <- IRanges::subsetByOverlaps(cigar.gr, lookup.region)
      cigar.grl <- split(cigar.gr, cigar.gr$qname)
      
      ## Construct plot ##
      ## Plot BAM alignements
      cigar.gr.df <- as.data.frame(cigar.gr)
      cigar.gr.df$level <- rep(1:length(unique(cigar.gr$qname)), lengths(cigar.grl))
      level.offset <- max(cigar.gr.df$level) + 1.5
      plt <- ggplot() + geom_rect(data=cigar.gr.df, aes(xmin=start, xmax=end, ymin=level, ymax=level+1, fill=cigar)) +
        scale_fill_manual(values = brewer.pal(n = 9, name = 'Set1'), drop=FALSE) +
        theme_bw() +
        ylab("") +
        xlab("Genomic Position (bp)") +
        ggtitle(regions.ID)
      ## Annotate read names
      text.label <- cigar.gr.df %>% group_by(qname) %>% summarise(x=min(start))
      text.label$y <- unique(cigar.gr.df$level)
      plt <- plt + geom_text(data=text.label, aes(x=x, y=y, label=qname), hjust=0, vjust=-0.5, inherit.aes = FALSE)
      ## Plot min & max range as gray bar
      max.region <- data.frame(start=min(c(start(cigar.gr), start(invDup.region))), 
                               end=max(c(end(cigar.gr), end(invDup.region))), y=level.offset)
      plt <- plt + geom_rect(data=max.region, aes(xmin=start, xmax=end, ymin=y+0.25, ymax=y+0.75), fill='gray', inherit.aes = FALSE)
      ## Plot region of interest as black bar
      invDup.region.df <- as.data.frame(invDup.region)
      invDup.region.df$y <- level.offset
      plt <- plt + geom_rect(data=invDup.region.df, aes(xmin=start, xmax=end, ymin=y, ymax=y+1), fill='black', inherit.aes = FALSE)
      plots[[length(plots) + 1]] <- plt
    } else {
      message("    No fragments assembled for this region!!!")
    }
  }
  if (!is.null(file)) {
    ## Export final plots in PDF
    message("Printing to PDF ...")
    grDevices::pdf(file, width=10, height=4)
    bquiet = lapply(plots, print)
    d <- grDevices::dev.off()
  }
  return(plots)
  message("DONE!!!")
}


#' Prepare NJ tree based on data matrix
#'
#' @param data.matrix ...
#' @param boot.iter ...
#' @return A \code{ggplot} object.
#' @importFrom ape nj boot.phylo
#' @author David Porubsky
#' @export
#'
plotDistanceTree <- function(data.matrix, boot.iter=10000) {
  ## Construct a tree
  tree <- ape::nj(X = dist(data.matrix))
  if (boot.iter > 0) {
    boot <- ape::boot.phylo(tree, data.matrix, function(x) nj(dist(x)), B = boot.iter)
    boot <- (boot/boot.iter)*100
    tree$node.label <- boot
  }
  offset <- max(tree$edge.length) + 100
  ## Plot phylogenetic tree
  plt <- ggplot(tree) + 
    geom_tree() + 
    theme_tree2() + 
    geom_tiplab() + 
    geom_nodelab(hjust = -0.1) +
    ggplot2::xlim(0, offset)
  return(plt)
}
