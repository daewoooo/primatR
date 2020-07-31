#' Function to calculate percentage overlap between query and the subject ranges.
#'
#' @param query.gr A \code{\link{GRanges-class}} object.
#' @param subject.gr A \code{\link{GRanges-class}} object. 
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export

getSegDupOverlaps <- function(query.gr, subject.gr) {
  query.gr$SDTrackPercOverlap <- 0
  query.gr$TotalUniqueBases <- 0
  query.gr$LongestUniqueBases <- 0
  hits <- GenomicRanges::findOverlaps(query.gr, subject.gr)
  
  for (i in 1:length(query.gr)) {
    query <- query.gr[i]
    subject <- subject.gr[subjectHits(hits)[queryHits(hits) == i]]
    if (length(subject) > 0) {
      ## Calculate % overlap with subject ranges
      overlap <- IRanges::intersect(query, subject)
      percOverlap <- (sum(width(overlap))/width(query))*100
      query.gr$SDTrackPercOverlap[i] <- percOverlap
      
      ## Get length of a range not overlaping with subject ranges
      nonoverlap <- IRanges::disjoin(c(query[,0], subject[,0]), ignore.strand=TRUE)
      query.unique <- subsetByOverlaps(nonoverlap, subject, invert=TRUE)
      if (length(query.unique) > 0) {
        query.gr$TotalUniqueBases[i] <- sum(width(query.unique))
        query.gr$LongestUniqueBases[i] <- max(width(query.unique))
      } 
    } else {
      query.gr$TotalUniqueBases[i] <- width(query)
      query.gr$LongestUniqueBases[i] <- width(query)
    }
  }
  return(query.gr)
}

#' Function to calculate percentage overlap between query and the subject ranges.
#'
#' @param query.gr A \code{\link{GRanges-class}} object.
#' @param subject.gr A \code{\link{GRanges-class}} object.
#' @param index An user defiend ID to disntiguish reported data columns.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export

getRangesOverlaps <- function(query.gr=NULL, subject.gr=NULL, index=NULL) {
  ## Define columns to report
  query.gr$PercOverlap <- 0
  query.gr$TotalUniqueBases <- 0
  query.gr$LongestUniqueBases <- 0
  
  ## Make sure subject.gr is non-overlapping
  subject.gr <- reduce(subject.gr)
  
  ## Find overlaps between query and subject ranges
  hits <- GenomicRanges::findOverlaps(query.gr, subject.gr)
  
  for (i in 1:length(query.gr)) {
    query <- query.gr[i]
    subject <- subject.gr[subjectHits(hits)[queryHits(hits) == i]]
    if (length(subject) > 0) {
      ## Calculate % overlap with subject ranges
      overlap <- IRanges::intersect(query, subject)
      percOverlap <- (sum(width(overlap))/width(query))*100
      query.gr$PercOverlap[i] <- percOverlap
      
      ## Get length of a range not overlaping with subject ranges
      nonoverlap <- IRanges::disjoin(c(query[,0], subject[,0]), ignore.strand=TRUE)
      query.unique <- subsetByOverlaps(nonoverlap, subject, invert=TRUE)
      if (length(query.unique) > 0) {
        query.gr$TotalUniqueBases[i] <- sum(width(query.unique))
        query.gr$LongestUniqueBases[i] <- max(width(query.unique))
      } 
    } else {
      query.gr$TotalUniqueBases[i] <- width(query)
      query.gr$LongestUniqueBases[i] <- width(query)
    }
  }
  
  ## Add index to newly added metacolumns (perc.overlap, toGR, idx)
  if (is.character(index) & nchar(index) > 0) {
    ## Add index to query ranges
    col.names <- c('PercOverlap','TotalUniqueBases','LongestUniqueBases')
    col.names.query <- names(mcols(query.gr))
    to.modif <- which(col.names.query %in% col.names)
    col.names.query[to.modif] <- paste0(col.names.query[to.modif], "_", index)
    names(mcols(query.gr)) <- col.names.query
  }
  return(query.gr)
}
