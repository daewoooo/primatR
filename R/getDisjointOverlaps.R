#' Function to group overlapping ranges based on percentage of the overlap.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param percTh A percentage threshold for a required overlap.
#' @return A \code{\link{GRanges-class}} object with extra meta-columns (idx, perc.overlap & group).
#' @author David Porubsky
#' @export
#' 
getDisjointOverlaps <- function(gr, percTh = 50) {
  
  ## Helper function
  returnMaxOverlap <- function(gr) {
    return( gr[which.max(gr$perc.overlap)] )
  }
  
  ## Initialize exported metadata columns
  gr$idx <- 1:length(gr) #unique id for each inputed range
  gr$perc.overlap <- 0
  gr$group <- 0 #group for overlapping ranges
  gr$sub.group <- 0 
  ## Repeat this for ranges with group assigned a zero value
  while (any(gr$group == 0)) {
    ## Select ranges to assign a non-zero group
    process.gr <- gr[gr$group == 0]
    ## Disjoin overlapping ranges into a set non-overlapping ranges
    gr.disjoin <- GenomicRanges::disjoin(process.gr)
    ## Find overlaps between original ranges and a set of non-overlapping ranges
    hits <- findOverlaps(process.gr, gr.disjoin)
    ## Calculate % overlaps between original and a set of non-overlapping ranges
    query.width <- width(process.gr[queryHits(hits)]) #Set of original(input) ranges
    subject.width <- width(gr.disjoin[subjectHits(hits)]) #Non-overlapping ranges corresponding to the original ranges
    perc.overlap <- (subject.width / query.width) * 100
    ## Prepare ranges for % overlap filtering
    gr2filt <- process.gr[queryHits(hits)]
    gr2filt$perc.overlap <- perc.overlap
    gr2filt$group <- subjectHits(hits) + max(gr$group)
    ## Filter ranges based on user defined percTh
    #gr.filt <- gr2filt[perc.overlap >= percTh]
    gr2filt.grl <- GenomicRanges::split(gr2filt, gr2filt$idx)
    gr2filt.grl <- endoapply(gr2filt.grl, returnMaxOverlap)
    gr.filt <- unlist(gr2filt.grl, use.names = FALSE)
    ## Split ranges passing percTh and recalibrate % of overlap
    gr.filt.grl <- GenomicRanges::split(gr.filt, gr.filt$group)
    gr.filt.grl <- endoapply(gr.filt.grl, recalcPercOverlap)
    new.gr <- unlist(gr.filt.grl, use.names = FALSE)
    ## Set ranges with % overlap less then percTh to zero
    #new.gr$perc.overlap[new.gr$perc.overlap < percTh] <- 0
    ## Set subgroup based on required percTh
    new.gr$sub.group <- paste0(new.gr$group,".", 1)
    new.gr$sub.group[which(new.gr$perc.overlap < percTh)] <- paste0(new.gr$group[which(new.gr$perc.overlap < percTh)],".", 2)
    ## Assign processed ranges to initial ranges
    gr[match(new.gr$idx, gr$idx)] <- new.gr[new.gr$idx %in% gr$idx]
  }
  return(gr)
}

#' Function to group overlapping ranges based on percentage of the overlap.
#' 
#' This function weights overlaps based on percentage of the overlap and the amount of overlapping ranges
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param percTh A percentage threshold for a required overlap.
#' @return A \code{\link{GRanges-class}} object with extra meta-columns (idx, perc.overlap & group).
#' @author David Porubsky
#' @export
#' 
getDisjointOverlapsWeighted <- function(gr, percTh = 50) {
  
  ## Helper function
  returnMaxWeigth <- function(gr) {
    return( gr[which.max(gr$weigth)] )
  }
  
  ## Initialize exported metadata columns
  gr$idx <- 1:length(gr) #unique id for each inputed range
  gr$perc.overlap <- 0
  gr$weigth <- 0 #overlap weight given percentage overlap and # of overlapping ranges
  gr$group <- 0 #group for overlapping ranges
  gr$sub.group <- 0 
  ## Repeat this for ranges with group assigned a zero value
  while (any(gr$group == 0)) {
    ## Select ranges to assign a non-zero group
    process.gr <- gr[gr$group == 0]
    ## Disjoin overlapping ranges into a set non-overlapping ranges
    gr.disjoin <- GenomicRanges::disjoin(process.gr)
    ## Find overlaps between original ranges and a set of non-overlapping ranges
    hits <- findOverlaps(process.gr, gr.disjoin)
    ## Remove hits overlapping only with one range
    disjoint.idx <- subjectHits(hits)
    non.unique.hits <- disjoint.idx[duplicated(disjoint.idx) | duplicated(disjoint.idx, fromLast=TRUE)]
    hits <- hits[subjectHits(hits) %in% non.unique.hits,]
    ## Calculate % overlaps between original and a set of non-overlapping ranges
    query.width <- width(process.gr[queryHits(hits)]) #Set of original(input) ranges
    subject.width <- width(gr.disjoin[subjectHits(hits)]) #Non-overlapping ranges corresponding to the original ranges
    perc.overlap <- (subject.width / query.width) * 100
    ## Remove overlaps below percTh
    #mask <- perc.overlap >= percTh
    #perc.overlap <- perc.overlap[mask]
    #hits <- hits[mask,]
    
    if (length(hits) > 0) {
      ## Prepare ranges for % overlap filtering
      gr2filt <- process.gr[queryHits(hits)]
      gr2filt$perc.overlap <- perc.overlap
      gr2filt$group <- subjectHits(hits) + max(gr$group)
      ## Calcualte weigth of each overlap based on % overlap and the number of overlapping ranges
      perc.overlap.perGroup <- split(gr2filt$perc.overlap, gr2filt$group)
      group.weigth <- sapply(perc.overlap.perGroup, function(x) sum(x)*length(x))
      gr2filt$weigth <- group.weigth[match(gr2filt$group, names(group.weigth))]
      ## Filter ranges based on user defined percTh
      #gr.filt <- gr2filt[perc.overlap >= percTh]
      gr2filt.grl <- GenomicRanges::split(gr2filt, gr2filt$idx)
      #gr2filt.grl <- endoapply(gr2filt.grl, returnMaxOverlap)
      gr2filt.grl <- endoapply(gr2filt.grl, returnMaxWeigth)
      gr.filt <- unlist(gr2filt.grl, use.names = FALSE)
      ## Split ranges passing percTh and recalibrate % of overlap
      gr.filt.grl <- GenomicRanges::split(gr.filt, gr.filt$group)
      gr.filt.grl <- endoapply(gr.filt.grl, recalcPercOverlap)
      new.gr <- unlist(gr.filt.grl, use.names = FALSE)
      ## Set groups with more than one range in a group to zero for overlap recalculation
      to.recalculate <- new.gr[new.gr$perc.overlap < percTh]
      if (length(to.recalculate) > 1 & length(new.gr) > length(to.recalculate)) {
        group.to.recalculate <- as.numeric(names(which(table(to.recalculate$group) > 1)))
        if (length(group.to.recalculate) > 0) {
          idx.to.recalculate <- to.recalculate$idx[to.recalculate$group %in% group.to.recalculate]
          new.gr <- new.gr[!new.gr$idx %in% idx.to.recalculate]
        }  
      }  
      ## Set ranges with % overlap less then percTh to zero
      #new.gr$perc.overlap[new.gr$perc.overlap < percTh] <- 0
      ## Set subgroup based on required percTh
      new.gr$sub.group <- paste0(new.gr$group,".", 1)
      if (any(new.gr$perc.overlap < percTh)) {
        n <- length(new.gr[new.gr$perc.overlap < percTh])
        sub.group.idx <- seq(from = 2, length.out = n)
        new.gr$sub.group[which(new.gr$perc.overlap < percTh)] <- paste0(new.gr$group[which(new.gr$perc.overlap < percTh)],".",  sub.group.idx)
      }
    } else {
      new.gr <- process.gr
      new.gr$group <- seq(max(gr$group) + 1, (max(gr$group) + length(new.gr)))
      new.gr$sub.group <- paste0(new.gr$group,".", 1)
    } 
    ## Assign processed ranges to initial ranges
    gr[match(new.gr$idx, gr$idx)] <- new.gr[new.gr$idx %in% gr$idx]
  }
  return(gr)
}


#' Function to calculate reciprocal overlap for a set of genomic ranges.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky

recalcPercOverlap <- function(gr) {
  if (length(gr) > 1) {
    disjoin.gr <- disjoin(gr)
    hits <- findOverlaps(disjoin.gr, gr)
    common.subrange.idx <- names( which(table(queryHits(hits)) == length(gr)) )
    gr$perc.overlap <- ( width(disjoin.gr[as.numeric(common.subrange.idx)]) / width(gr) ) * 100
  } else {
    gr$perc.overlap <- 0
  }  
  return(gr)
}

