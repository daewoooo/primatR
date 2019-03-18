#' Filter and plot significant connections between read pairs
#' 
#' This function takes condidate links and filter all spurious links. Also this function exports table of links 
#' as well as plot links onto the genome reference.
#'
#' @param gr.links A \code{\link{GRanges-class}} object with stored links between read pairs.
#' @inheritParams readPairsAsLinks
#' @author David Porubsky
#' @export

processReadLinks <- function(gr.links, min.reads = 10, chromosomes = NULL, bsgenome = NULL) {
  
  ## Helper functions
  getUniqueLinks <- function(readIDs) {
    comb.read.IDs <- t(combn(length(readIDs), 2))
    counts <- list()
    for (i in 1:nrow(comb.read.IDs)) { 
      pair <- comb.read.IDs[i,]
      count <-  length(intersect(readIDs[[pair[1]]], readIDs[[pair[2]]]))
      counts[[i]] <- c(pair, count)
    }
    counts <- do.call(rbind, counts)
    return(counts)
  }

  ## Get number of assessed breaks from the link object
  assessed.breaks <- gr.links$assessed.breaks
  ## Get list of procesed bam files
  assessed.bams <- gr.links$assessed.bams
  
  ## Process links
  links.gr <- gr.links$links
  links.gr$break.chr <- sapply(links.gr$ID, function(x) gsub(strsplit(x, "-")[[1]][1], pattern = 'chr', replacement = ''))
  
  ## Add x_coord by selecting larger number being the start of the link
  links.gr$x <- pmax(start(links.gr), start(links.gr$to.gr))
  links.gr$xend <- pmin(start(links.gr), start(links.gr$to.gr))
  ## Remove links with the same start and end position
  links.gr$link <- paste0(links.gr$x, "_", links.gr$xend)
  links.gr <- links.gr[!duplicated(links.gr$link),]
  links.gr <- links.gr[links.gr$x != links.gr$xend,]
  ## Filter only user defined chromosomes
  mask <- as.vector(seqnames(links.gr)) %in% chromosomes & as.vector(seqnames(links.gr$to.gr)) %in% chromosomes
  links.gr <- links.gr[mask,]

  ## Impose chromosome levels to observed data
  links.gr$y <- as.character(seqnames(links.gr))
  links.gr$yend <- as.character(seqnames(links.gr$to.gr))
  links.gr$y <- match(links.gr$y, as.character(chromosomes))
  links.gr$yend <- match(links.gr$yend, as.character(chromosomes))
  ## Make sure that at least one break resides on chromosome of origin
  links.gr <- links.gr[links.gr$break.chr == links.gr$y | links.gr$break.chr == links.gr$yend]
  
  ## Add links between chromosomes
  chr.from <- pmax(links.gr$y, links.gr$yend)
  chr.to <- pmin(links.gr$y, links.gr$yend)
  links.gr$chr.link <- paste0(chr.from, "_", chr.to)
  links.gr$chr.link.perBreak <- paste0(links.gr$ID, "_", links.gr$chr.link)
  
  ## Get number of unique reads per link
  qnames.list <- split(links.gr$qname, links.gr$chr.link.perBreak)
  qnames.counts.perLink <- sapply(qnames.list, countUniqueReadIDs)
  ## Filter links with minimal unique read support
  signif.links <- names(qnames.counts.perLink)[qnames.counts.perLink >= min.reads]
  links.gr.filt <- links.gr[links.gr$chr.link.perBreak %in% signif.links,]
  
  ## Get number of intra- versus inter-chromosomal links
  ## Export siginificant links
  total.breaks <- length(unique(links.gr.filt$ID))
  breaks.grl <- split(links.gr.filt, links.gr.filt$ID)
  intra.links.count <- 0
  inter.links.count <- 0
  intra.links.list <- GRangesList()
  inter.links.list <- GRangesList()
  final.links <- GRangesList()
  for (i in seq_along(breaks.grl)) {
    #print(i)
    break.gr <- breaks.grl[[i]]
    roi.ID <- strsplit(break.gr$qname[1], "__")[[1]][3]
    links.grl <- split(break.gr, break.gr$chr.link)

    for (j in seq_along(links.grl)) {
      #link.df <- NULL
      link.processed.gr <- NULL
      link.gr <- links.grl[[j]]
      ## Get collapsed genomic regions of selected links
      link.ID <- unique(link.gr$ID)
      all.ranges <- c(link.gr[,'qname'], link.gr$to.gr[,'qname'])
      strand(all.ranges) <- "*"
      regions <- reduce(all.ranges, min.gapwidth=10000) #Merge regions with minimal gap 10kb
      ## Remove regions with low link count
      overlaps <- countOverlaps(regions, all.ranges)
      regions <- regions[overlaps >= min.reads]
      ## Filter out reads from regions with low link count
      all.ranges <- subsetByOverlaps(all.ranges, regions)
      mask1 <- findOverlaps(link.gr, regions)
      mask2 <- findOverlaps(link.gr$to.gr, regions)
      mask <- intersect(queryHits(mask1), queryHits(mask2))
      link.gr <- link.gr[mask]
      ## Store all links in a list
      final.links[[length(final.links) + 1]] <- link.gr
      
      ## Find links with set minimal read support (min.reads)
      if (length(regions) > 1) {
        hits <- findOverlaps(all.ranges, regions)
        read.names <- split(all.ranges$qname, subjectHits(hits))
        link.counts <- getUniqueLinks(readIDs = read.names)
        max.link.idx <- which.max(link.counts[,3])  ## TODO check within chrom links being exported!!!
        #if (link.counts[max.link.idx, 3] >= min.reads) {
        region.IDs <- link.counts[max.link.idx, c(1:2)]
        link.count <- link.counts[max.link.idx, 3]
        regions.new <- regions[region.IDs]
        
        ## Store processed links
        link.processed.gr <- regions.new[1]
        link.processed.gr$to.gr <- regions.new[2]
        link.processed.gr$roi.ID <- roi.ID
        link.processed.gr$break.ID <- link.ID
        link.processed.gr$x <- link.gr$y[1]
        link.processed.gr$xend <- link.gr$yend[1]
        link.processed.gr$link.count <- link.count
      } else {
        link.count <- countOverlaps(regions, all.ranges)
        ## Store processed links
        link.processed.gr <- regions
        link.processed.gr$to.gr <- regions
        link.processed.gr$roi.ID <- roi.ID
        link.processed.gr$break.ID <- link.ID
        link.processed.gr$x <- link.gr$y[1]
        link.processed.gr$xend <- link.gr$yend[1]
        link.processed.gr$link.count <- link.count
      } 
          
      ## Store inter- and intra-links in separate objects
      if (!is.null(link.processed.gr) & length(link.gr) > 0) {
        #print(link.df)
        chroms <- strsplit(unique(link.gr$chr.link), "_")[[1]]
        if (chroms[1] == chroms[2]) {
          intra.links.count <- intra.links.count + 1
          link.gr$valid.ID <- paste0('valid', length(intra.links.list) + 1)
          intra.links.list[[length(intra.links.list) + 1]] <-  link.processed.gr
        } else {
          inter.links.count <- inter.links.count + 1
          link.gr$valid.ID <- paste0('valid', length(inter.links.list) + 1)
          inter.links.list[[length(inter.links.list) + 1]] <-  link.processed.gr
        }
      }  
    }
  }  

  ## Export final data objects
  #intra.links <- do.call(c, intra.links.list)
  #inter.links <- do.call(c, inter.links.list)
  intra.links <- unlist(intra.links.list)
  inter.links <- unlist(inter.links.list)
  final.links <- unlist(final.links, use.names = FALSE)
  missed.bams <- assessed.bams[!assessed.bams %in% unique(final.links$bam.name)]
  
  ## Collapse overlapping regions
  intra.links.hits <- findOverlaps(intra.links, drop.self=TRUE, select = 'first')
  inter.links.hits <- findOverlaps(inter.links, drop.self=TRUE, select = 'first')
  
  return(list(links.gr=final.links, intra.links=intra.links, inter.links=inter.links, assessed.breaks=assessed.breaks, missed.bams=missed.bams))
}  
    

#' Plot significant connections between read pairs
#' 
#' This function takes condidate filtered links exported from function 'processReadLinks'
#'
#' @param links A \code{\link{GRanges-class}} object with stored links between read pairs.
#' @inheritParams readPairsAsLinks
#' @importFrom scales comma
#' @author David Porubsky
#' @export
plotLinks <- function(links=NULL, chromosomes=NULL, index=NULL, bsgenome=NULL) {
  links.intra <- links$intra.links
  links.inter <- links$inter.links
  
  ## Create non-redundant links
  # inter.hits1 <- findOverlaps(links.inter, drop.self=TRUE, select = "first")
  # inter.hits2 <- findOverlaps(links.inter$to.gr, drop.self=TRUE, select = "first")
  # mask.inter <- inter.hits1 == inter.hits2
  # mask.inter[is.na(mask.inter)] <- FALSE
  # intra.hits1 <- findOverlaps(links.intra, drop.self=TRUE, select = "first")
  # intra.hits2 <- findOverlaps(links.intra$to.gr, drop.self=TRUE, select = "first")
  # mask.intra <- intra.hits1 == intra.hits2
  # mask.intra[is.na(mask.intra)] <- FALSE
  
  links.inter.df <- as.data.frame(links.inter)
  links.intra.df <- as.data.frame(links.intra)
  
  ## Prepare genome-wide ideogram
  seq.len <- seqlengths(bsgenome)[chromosomes]
  ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
  ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chromosomes)
  ideo.df$levels <- 1:length(seq.len)
  ideo <- ggplot(ideo.df) + geom_linerange(aes(x=levels, ymin=0, ymax=length), size=1, color='black')
  
  ## Add inter chromosomal links to the ideogram
  plt <- ideo +
    geom_curve(data=links.inter.df, aes(x = x, y = start, xend = xend, yend = to.gr.end), color="red", curvature = 0.25) +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(breaks = ideo.df$levels, labels = ideo.df$seqnames) +
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    theme_bw() +
    ylab("Genomic position (bp)") +
    xlab("Chromosome")
  
  ## Add inter chromosomal links to the ideogram
  plt <- plt +
    geom_curve(data=links.intra.df, aes(x = x, y = start, xend = xend, yend = to.gr.end), color="orange", curvature = 0.25)
  
  ## Add points of all intra-chromosomal links
  intra.points.df <- data.frame(x=c(links.intra.df$x, links.intra.df$xend), y=c(links.intra.df$start, links.intra.df$to.gr.end), roi.ID=rep(links.intra.df$roi.ID, 2))
  plt <- plt +  
    geom_point(data=intra.points.df, aes(x = x, y = y), color="orange", shape=95, size=5, inherit.aes = FALSE)
  
  ## Add some count statistics
  all.assessed.breaks <- links$assessed.breaks
  intra.links.count <- length(links.intra)
  inter.links.count <- length(links.inter)
  #unique.inter <- length(unique(links.inter$break.ID))
  #unique.intra <- length(unique(links.intra$break.ID))
  inter.and.intra.breaks <- length(intersect(links.intra$break.ID, links.inter$break.ID))
  only.intra <- length(unique(setdiff(links.intra$break.ID, links.inter$break.ID)))
  only.inter <- length(unique(setdiff(links.inter$break.ID, links.intra$break.ID)))
  multi.inter <- length(links.inter[duplicated(links.inter$break.ID)])
  
  plt <- plt + annotate(geom="text", x=23, y=250000000, label=paste0("All assessed breaks: ", all.assessed.breaks, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-10000000, label=paste0("All intra-chromosomal breaks: ", intra.links.count, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-20000000, label=paste0("All inter-chromosomal breaks: ", inter.links.count, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-30000000, label=paste0("Only intra-chromosomal breaks: ", only.intra, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-40000000, label=paste0("Only inter-chromosomal breaks: ", only.inter, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-50000000, label=paste0("Both intra & inter breaks: ", inter.and.intra.breaks, "\n"), color="black", hjust=1)
  plt <- plt + annotate(geom="text", x=23, y=250000000-60000000, label=paste0("Multi-location inter-chromosomal breaks: ", multi.inter, "\n"), color="black", hjust=1)
  
  ## Add plot title
  if (!is.null(index) & nchar(index) > 0) {
    plt <- plt + ggtitle(index)
  }
  
  return(plt)
}