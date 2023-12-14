#' This function will take in a \code{\link{GRanges-class}} object of genomic ranges and enumerate the number of overlapping bases
#' with other set of user defined 'query' genomic ranges.
#'
#' @param ranges A \code{\link{GRanges-class}} object of genomic regions to get the number of overlapping bases with 'query.ranges'.
#' @param query.ranges A \code{\link{GRanges-class}} object of genomic regions for which one want to report overlaps. 
#' @param index A user defined name of a column where the number of overlapping bases between 'ranges' and 'query.ranges' will be reported.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
reportOverlapBases <- function(ranges=NULL, query.ranges=NULL, index=NULL) {
  if (!is.null(query.ranges) & !is.null(ranges)) {
    if (class(ranges) == "GRanges" & class(query.ranges) == "GRanges") {
      ## Make sure strand is not considered
      strand(ranges) <- '*'
      strand(query.ranges) <- '*'
      ## Get disjoint ranges between query and user.defined set of ranges
      disj.gr <- suppressWarnings( GenomicRanges::disjoin(c(ranges[,0], query.ranges[,0])) )
      ## Get disjoint ranges that are overlapping with query.ranges
      disj.query.gr <- IRanges::subsetByOverlaps(disj.gr, query.ranges)
      ## Get disjoint ranges overlapping with ranges of interest
      disj.query.roi.gr <- IRanges::subsetByOverlaps(disj.query.gr, ranges)
      ## Split by regions of interest
      hits <- IRanges::findOverlaps(ranges, disj.query.roi.gr)
      disj.query.roi.grl <- split(disj.query.roi.gr[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))
      query.bases <- sapply(disj.query.roi.grl, function(gr) sum(width(reduce(gr))))
      ## Add overlapping bases counts
      if (!is.null(index) & nchar(index) > 0) {
        new.col.idx <- ncol(GenomicRanges::mcols(ranges)) + 1
        GenomicRanges::mcols(ranges)[new.col.idx] <- 0
        colnames(GenomicRanges::mcols(ranges))[new.col.idx] <- index
        GenomicRanges::mcols(ranges)[new.col.idx][unique(S4Vectors::queryHits(hits)),] <- query.bases
      } else {
        ranges$query.bases <- 0
        ranges$query.bases[unique(S4Vectors::queryHits(hits))] <- query.bases
      } 
    }
  }
  return(ranges)
}
