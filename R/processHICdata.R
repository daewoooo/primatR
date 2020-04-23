#' Construct and plot HIC contact matrix.
#'
#' This function imports HIC reads and plot contact matrix based on user defined resolution. 
#' 
#' @param resolution A \code{vector} of genomic bin sizes for which contact counts will be constructed and plotted.
#' @param highlight.pos A \code{vector} of genomic positions to be highlighted as horizontal and vertical lines.
#' @param title A \code{character} string to be used as a title for each plot.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @inheritParams bam2GRanges
#' @inheritParams plotNucmerCoords
#' @author David Porubsky
#' @export
plotHICregional <- function(bamfile, region = NULL, min.mapq = 10, resolution = 50000, highlight.pos = NULL, bsgenome = NULL, title=NULL) {
  
  ## Read in paired-end reads
  message("Loading reads for region: ", as.character(region), " ...")
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile,
                                                                   index=paste0(bamfile,'.bai'),
                                                                   param=Rsamtools::ScanBamParam(tag="XA",
                                                                                                 which=region,
                                                                                                 what=c('mapq','flag','mrnm','mpos'), 
                                                                                                 mapqFilter = min.mapq,
                                                                                                 flag=Rsamtools::scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE)
                                                                   ), 
                                                                   use.names = TRUE
  )
  )
  ## Convert reads into read pairs
  data.raw <- GenomicAlignments::makeGAlignmentPairs(x = data.raw)
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
  # if (!is.null(min.mapq)) {
  #   mask <- data.first$mapq >= min.mapq & data.last$mapq >= min.mapq
  #   data.first <- data.first[mask]
  #   data.last <- data.last[mask]
  # }
  
  raw.data <- list()
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
    link.df <- data.frame('M1'=as.data.frame(bins[pairs[,1]])[,1:3], 'M2'=as.data.frame(bins[pairs[,2]])[,1:3], value=link.counts, idx1 = pairs[,1], idx2 = pairs[,2])
    ## Store raw data
    raw.data[[length(raw.data) + 1]] <- link.df
    
    ## Plot binned contacts
    plt <- ggplot(link.df[link.df$value > 0,]) + 
      geom_rect(aes(xmin=M1.start, xmax=M1.end, ymin=M2.start, ymax=M2.end, fill=value))
    if (!is.null(highlight.pos)) {
      highlight.pos <- sort(highlight.pos)
      plt <- plt + geom_vline(xintercept = highlight.pos[1], color='red') +
        geom_hline(yintercept = highlight.pos[2], color='red')
    }
    ## Add title if defined to the resolution label
    if (!is.null(title) & nchar(title > 0)) {
      plt.label <- paste0(title, "\nResolution: ", res, "bp")
    } else {
      plt.label <- paste0("Resolution: ", res, "bp")
    }
    ## Construct final plot
    plt <- plt + scale_fill_gradient(low = "white", high = "black", trans = 'log', name="Links (log)") +
      scale_x_continuous(labels = comma) +
      scale_y_continuous(labels = comma) +
      xlab("Genomic position (bp)") +
      ylab("Genomic position (bp)") +
      annotate(geom="text", x=max(link.df$M1.start), y=min(link.df$M1.start), label=plt.label, hjust=1, vjust=-0.5) +
      theme(aspect.ratio=1) + #Make sure plot is always a square
      theme_bw()
    plots[[length(plots) + 1]] <- plt 
  }
  #message("DONE!!!")
  ## Return plots and raw data used for plotting
  return(list(plots = plots, raw.data = raw.data))
}  



#' Plot horizontal HIC contact matrix.
#'
#' This function imports HIC contact values and plots horizontal map of HIC contacts. 
#' 
#' @param HICcontacts.df A \code{data.frame} with the following columns: M1.seqnames, M1.start, M1.end, M2.seqnames, M2.start, M2.end, value, idx1, idx2
#' @param yHeightRatio A size ratio between x and y axis. Default: 1/3, what means that the top third of y-axis values will be removed.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' @author David Porubsky
#' @export
plotHICcontactTable <- function(HICcontacts=NULL, yHeightRatio=1/3, sd.track=NULL, user.track=NULL, gene.track=NULL, title=NULL) {  
  ## Get region boundaries
  xstart <- min(HICcontacts$M1.start)
  xend <- max(HICcontacts$M2.end)
  region.width <- xend - xstart
  #binNum <- (xend - xstart)/binSize
  ## Add unique bin id to each 
  HICcontacts$binID <- paste0('bin', 1:nrow(HICcontacts))
  ## Get differences between matrix indices
  bin.diffs <- HICcontacts$idx2 - HICcontacts$idx1
  ## Get binsize of used genomic bins
  bin.size <- (HICcontacts$M1.end[1] - HICcontacts$M1.start[1]) + 1
  ## Get position of bin centers on x-axis
  centers <- rowMeans(
    data.frame(
      HICcontacts$M1.start,
      HICcontacts$M2.end ) )
  ## Get bin positions on y-axis
  maxys <- 1 + (bin.diffs - 1) * 0.5
  minys <- maxys - 1
  ## Construct data.frame for plotting
  coords.df <- rbind(
    data.frame(
      x=centers - (bin.size/2),
      y=rowMeans( data.frame(maxys, minys) ),
      value=HICcontacts$value,
      bin=HICcontacts$binID),
    data.frame(
      x=centers,
      y=maxys,
      value=HICcontacts$value,
      bin=HICcontacts$binID),
    data.frame(
      x=centers + (bin.size/2),
      y=rowMeans( data.frame(maxys, minys) ),
      value=HICcontacts$value,
      bin=HICcontacts$binID),
    data.frame(
      x=centers,
      y=minys,
      value=HICcontacts$value,
      bin=HICcontacts$binID) )
  
  ## Get center of each bin on x-axis
  coords.df$x <- coords.df$x - bin.size/2
  ## Blunt ends of x-axis
  xstart <- xstart - bin.size/2
  xend <- xend - bin.size/2
  ## Scale y-axis to bin.size and set limit based on user defined retion of x and y size
  coords.df$y <- pmax( 0, coords.df$y )
  coords.df$y <- coords.df$y * bin.size
  ylim2 <- region.width * yHeightRatio
  
  ## Set ploting theme
  my.theme <- theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())
  ## Make the plot
  plt <- ggplot(coords.df[coords.df$value > 0,], aes(x = x, y = y)) +
    geom_polygon(aes(fill = value, group = bin)) +
    scale_fill_gradient(low = "white", high = "black", trans = 'log', name="Links (log)") + 
    #coord_fixed( ratio=1, xlim=c(xstart, xend), ylim=c(0, ylim2) ) +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(labels = comma) +
    xlab("Genomic Position (bp)") +
    ylab("") +
    theme_bw() +
    my.theme
  
  offset <- -1
  ## Add SD annotation
  if (!is.null(sd.track) & length(sd.track) > 0) {
    sd.track <- subsetByOverlaps(sd.track, region)
    sd.track <- reduce(sd.track)
    sd.track.df <- as.data.frame(sd.track)
    sd.track.df$ID <- 'SDs'
    ylim2.minus <- offset
    plt <- plt + geom_segment(data=sd.track.df, aes(x = start, y = offset, xend = end, yend = offset, color=ID), size=1.5, inherit.aes = FALSE)
  }
  ## Add user annotation tracks if defined
  if (!is.null(user.track) & length(user.track) > 0) {
    user.track$level <- -disjointBins(user.track)
    user.track$level <- user.track$level + offset
    annot.df <- as.data.frame(user.track)
    scale <- ylim2 * 0.01
    ylim2.minus <- min(user.track$level * scale)
    #plt + new_scale(new_aes = 'fill') +
    #  geom_rect(data=annot.df, aes(xmin=start, xmax=end, ymin=level, ymax=level*scale, fill=gen), inherit.aes = FALSE)
    plt <- plt + geom_segment(data=annot.df, aes(x = start, y = level*scale, xend = end, yend = level*scale, color=gen), size=1.5, inherit.aes = FALSE)
    offset <- min(user.track$level)
  }
  ## Add gene annotation tracks if defined
  if (!is.null(gene.track) & length(gene.track) > 0) {
    gene.track$level <- -disjointBins(gene.track)
    gene.track$level <- gene.track$level + offset
    #annot.df <- as.data.frame(gene.track)
    tmp.df <- as.data.frame(gene.track)
    scale <- ylim2 * 0.03
    ylim2.minus <- min(gene.track$level * scale)
    plt <- plt + geom_segment(data=tmp.df, aes(x = start, y = level*scale, xend = end, yend = level*scale), size=1.5, inherit.aes = FALSE)
    plt <- plt + geom_text(data=tmp.df, aes(x = start, y = level*scale, label=gene.name), size=1.5, inherit.aes = FALSE, hjust=1, position = position_jitter())
  }
  ## Add y-axis limits and colors scale for user defined ranges
  plt <- plt + coord_fixed(ratio=1, xlim=c(xstart, xend), ylim=c(ylim2.minus, ylim2)) +
    scale_color_manual(values = c(HET='dodgerblue3', HOM='firebrick3', SDs='orange'), name="")
  
  ## Add title if defined to the resolution label
  if (!is.null(title) & nchar(title > 0)) {
    plt <- plt + ggtitle(toupper(title))
  }
  
  ## Return the plot
  return(plt)
}



#' Construct HIC contact matrix.
#'
#' This function imports HIC reads and construct contact matrix based on 
#' user defined resolution or a number of required genomic bins.
#' 
#' @param nbins A number of required bins to split selected genomic region into.
#' @param resolution A \code{vector} of genomic bin sizes for which contact counts will be constructed and plotted.
#' @return A \code{data.frame} object with HiC contact counts between two genomic region in each row.
#' @inheritParams bam2GRanges
#' @author David Porubsky
#' @export
getHIContactCounts <- function(bamfile, region = NULL, min.mapq = 10, nbins=NULL, resolution = NULL, bsgenome = NULL) {
  
  ## Read in paired-end reads
  message("Loading reads for region: ", as.character(region), " ...")
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile,
                                                                   index=paste0(bamfile,'.bai'),
                                                                   param=Rsamtools::ScanBamParam(tag="XA",
                                                                                                 which=region,
                                                                                                 what=c('mapq','flag','mrnm','mpos'), 
                                                                                                 mapqFilter = min.mapq,
                                                                                                 flag=Rsamtools::scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE)
                                                                   ), 
                                                                   use.names = TRUE
  )
  )
  ## Convert reads into read pairs
  data.raw <- GenomicAlignments::makeGAlignmentPairs(x = data.raw)
  ## Convert data into GRanges objects
  data.first <- as(GenomicAlignments::first(data.raw), 'GRanges') #mate1
  data.last <- as(GenomicAlignments::last(data.raw), 'GRanges') #mate2
  
  ## Filter reads that spans user-defined region
  hits.first <- findOverlaps(query = data.first, subject = region)
  hits.last <- findOverlaps(query = data.last, subject = region)
  mask <- intersect(queryHits(hits.first), queryHits(hits.last))
  data.first <- data.first[mask]
  data.last <- data.last[mask]
  
  ## Calculate binsize resolution
  if (nbins > 0 & is.null(resolution)) {
    app.resolution <- width(region) / nbins
    resolution <- plyr::round_any(app.resolution, 1000)
  } else if (nbins > 0 & resolution > 0) {
    message("    Warning: both 'nbins' and 'resolution' parameters defined, using 'nbins' to bin the genome!!!")
    app.resolution <- width(region) / nbins
    resolution <- plyr::round_any(app.resolution, 1000)
  } else {
    message("    Warning: please define 'nbins' or 'resolution' parameters to bin the genome!!!")
  }
  message("    Constructing contact matrix for binsize ", resolution, "bp ...")
  
  ## Make genomic bins
  bins <- makeBins(bsgenome = bsgenome, chromosomes = as.character(seqnames(region)), binsize = resolution, stepsize = resolution)
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
  link.df <- data.frame('M1'=as.data.frame(bins[pairs[,1]])[,1:3], 'M2'=as.data.frame(bins[pairs[,2]])[,1:3], value=link.counts, idx1 = pairs[,1], idx2 = pairs[,2])
  
  return(link.df)
} 