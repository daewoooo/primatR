#' Function to calculate percentage of reciprocal overlap between two sets of ranges.
#'
#' @param query A \code{\link{GRanges-class}} object.
#' @param subject A \code{\link{GRanges-class}} object.
#' @param thresh A percentage threshold for a required overlap.
#' @param report Set to 'both' if you want to return merged query and subject ranges or select only one of them. 
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky

getReciprocalOverlaps <- function(query, subject, thresh=50, report='both') {
  
  hits <- findOverlaps(query, subject)
  
  query$perc.overlap <- 0
  query$toGR <- query[,0]
  query$idx <- paste0(1:length(query), ".1")
  subject$perc.overlap <- 0
  subject$toGR <- subject[,0]
  subject$idx <- paste0(1:length(subject), ".2")
  for (i in seq_along(hits)) {
    ## Calculate % overlaps for a pair of ranges
    hit <- hits[i]
    percOverlap <- gr2grPercOverlap(gr1=query[queryHits(hit)], gr2=subject[subjectHits(hit)])
    
    ## Record only max overlap ranges
    if (query[queryHits(hit)]$perc.overlap < percOverlap) {
      query[queryHits(hit)]$perc.overlap <- percOverlap
      query[queryHits(hit)]$toGR <- subject[subjectHits(hit)][,0]
    }
    
    if (subject[subjectHits(hit)]$perc.overlap < percOverlap) {
      subject[subjectHits(hit)]$perc.overlap <- percOverlap
      subject[subjectHits(hit)]$toGR <- query[queryHits(hit)][,0]
    }  
    
    ## Assign the same index to ranges that meet reciprocal overlap threshold  
    if (percOverlap >= thresh) {
      query[queryHits(hit)]$idx <- queryHits(hit)
      subject[subjectHits(hit)]$idx <- queryHits(hit)
    }
  }
  
  if (report == 'query') {
    return(query)
  } else if (report == 'subject') {
    return(subject)
  } else {
    return(c(query, subject))
  } 
}


#' Function to calculate percentage overlap between two ranges (gr1 & gr2).
#'
#' @param gr1 A \code{\link{GRanges-class}} object.
#' @param gr2 A \code{\link{GRanges-class}} object. 
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky

gr2grPercOverlap <- function(gr1, gr2) {
  intersection <- IRanges::intersect(gr1, gr2)
  if (length(intersection) > 0) {
    gr1.overlap <- (width(intersection)/width(gr1))*100
    gr2.overlap <- (width(intersection)/width(gr2))*100
    recip.overlap <- min(gr1.overlap, gr2.overlap)
  } else {
    recip.overlap <- 0
  }  
  return(recip.overlap)
}