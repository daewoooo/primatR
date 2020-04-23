#' Function to subtract ranges that flank regions of interest
#'
#' @param gr A \code{\link{GRanges-class}} object of ranges from which remove.gr are removed.
#' @param remove.gr A \code{\link{GRanges-class}} object of ranges to subtract from gr.
#' @param mode Select 'flanks' to subtract only flanking regions or 'all' to subtract all regions what creates disjoined ranges.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
subtractRegions <- function(gr = NULL, remove.gr = NULL, mode = 'flanks') {
  ## Get set of non-overlapping ranges between gr of interest and mask - remove.gr
  disjoin.gr <- disjoin(c(gr[,0], remove.gr[,0]))
  ## Get only non-overlapping ranges that overlaps with gr of interest
  disjoin.gr <- subsetByOverlaps(disjoin.gr, gr)
  ## Remove non-overlapping ranges that overlaps with mask - remove.gr
  unique.gr <- subsetByOverlaps(disjoin.gr, remove.gr, invert = TRUE)
  ## Find overlaps between unique ranges and gr of interest
  hits <- findOverlaps(gr, unique.gr) 
  ## Merge splitted unique ranges that belong to the same gr of interest
  unique.gr.collapse <- unique.gr[subjectHits(hits)]
  unique.gr.collapse$ID <- queryHits(hits)
  ## Subtract only flanking regions overlapping with SDs
  if (mode == 'flanks') {
    ## Collapse ranges longer or equal 500 bp
    unique.gr.collapse <- unique.gr.collapse[width(unique.gr.collapse) >= 500]
    unique.merged.gr <- collapseBins(unique.gr.collapse, id.field = 1)
    return(unique.merged.gr)
  ## Report the single longest range non-overlapping with SDs 
  } else if (mode == 'longest') {
    unique.gr.collapse.grl <- split(unique.gr.collapse, unique.gr.collapse$ID)
    unique.gr.collapse.grl <- endoapply(unique.gr.collapse.grl, function(x) reduce(x)[which.max(width(reduce(x)))])
    return(unlist(unique.gr.collapse.grl, use.names = FALSE))
  } else {
    return(unique.gr.collapse)
  } 
}
