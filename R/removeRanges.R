#' Function to subtract set of genomic ranges from other genomic ranges
#'
#' @param gr A \code{\link{GRanges-class}} object of ranges from which remove.gr are removed.
#' @param remove.gr A \code{\link{GRanges-class}} object of ranges to subtract from gr.
#' @param collapse.gaps Size of the gaps between ranges to be merged.
#' @param min.region Minimum size regions to be reported.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
removeRanges <- function(gr=NULL, remove.gr=NULL, collapse.gaps=0, min.region=0) {
  ## Helper function
  returnMaxRegion <- function(gr, collapse.gaps=0, min.region=0) {
    gr <- keepSeqlevels(gr, unique(seqnames(gr)), pruning.mode = 'coarse')
    if (collapse.gaps > 0) {
      gr.gaps <- gaps(gr, start = min(start(gr)), end = max(end(gr)))
      gr.gaps <- gr.gaps[strand(gr.gaps) == '*']
      gr.gaps <- gr.gaps[width(gr.gaps) <= collapse.gaps]
      if (length(gr.gaps) > 0) {
        gr <- reduce(c(gr, gr.gaps))
        #gr <- gr[which.max(width(gr))]
      } else {
        gr <- reduce(gr)
        #gr <- gr[which.max(width(gr))]
      }
    } else {
      gr <- reduce(gr)
      #gr <- gr[which.max(width(gr))]
    }
    if (min.region > 0) {
      return(gr[width(gr) >= min.region])
    } else {
      return(gr)
    }  
  }
  
  ## Get set of non-overlapping ranges between gr of interest and mask - remove.gr
  disjoin.gr <- disjoin(c(gr[,0], remove.gr[,0]))
  ## Remove non-overlapping ranges that overlaps with mask - remove.gr
  unique.gr <- subsetByOverlaps(disjoin.gr, remove.gr, invert = TRUE)
  ## Get only non-overlapping ranges that overlaps with gr of interest
  hits <- findOverlaps(unique.gr, gr)
  unique.grl <- split(unique.gr[queryHits(hits)], subjectHits(hits))
  ## Report longest continous range that do not overlap with SD
  unique.reduced.grl <- endoapply(unique.grl, function(x) returnMaxRegion(x, collapse.gaps = collapse.gaps, min.region = min.region))
  unique.reduced.gr <- unlist(unique.reduced.grl, use.names = FALSE)
  return(unique.reduced.gr)
}
