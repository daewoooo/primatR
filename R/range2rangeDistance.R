#' Function to calculate distances to user defined ranges
#' 
#' This function calculates distance of each breakpoint of submitted \code{\link{GRanges-class}} object to its nearest user defined range 
#' defined in 'userTrack'.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param userTrack A \code{\link{GRanges-class}} object.
#' @param allow.overlap Set to \code{TRUE} if overlaps between ranges stored in 'gr' with ranges stored in 'userTrack' should be allowed
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
range2rangeDistance <- function(gr, userTrack, allow.overlap=FALSE) {
  ## Export start and end position as a separate GRanges objects
  gr.starts <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=start(gr), end=start(gr)+1))
  gr.ends <- GRanges(seqnames=seqnames(gr), ranges=IRanges(start=end(gr)-1, end=end(gr)))
  
  ## Reduce overlapping ranges
  strand(userTrack) <- "*" #Remove strand information
  userTrack.collapsed <- reduce(userTrack)
  
  ## Find closest range to left and the right breakpoint
  nearest.left.idx <-  follow(gr.starts, userTrack.collapsed) #get nearest range to the left from the INV start (upstream)
  nearest.right.idx <-  precede(gr.ends, userTrack.collapsed) #get nearest range to the right from the INV end (downstream)
  
  ## Resolve NA values
  left.na <- which(is.na(nearest.left.idx))
  right.na <- which(is.na(nearest.right.idx))
  ## Assign range 1 to all NA values
  nearest.left.idx[left.na] <- 1
  nearest.right.idx[right.na] <- 1
  ## Select nearest userTrack ranges
  nearest.left <- userTrack.collapsed[nearest.left.idx]
  nearest.right <- userTrack.collapsed[nearest.right.idx]
  ## Change NA values back to their own value
  nearest.left[left.na] <- gr.starts[left.na]
  nearest.right[right.na] <- gr.ends[right.na]
  ## Get distance to the nearest SD
  left.distances <- distance(nearest.left, gr.starts)
  right.distances <- distance(nearest.right, gr.ends)
  ## Set distances for NA values to be a high number
  left.distances[left.na] <- -1
  right.distances[right.na] <- -1
  
  if (allow.overlap) {
    left.hits <- findOverlaps(gr.starts, userTrack.collapsed)
    right.hits <- findOverlaps(gr.ends, userTrack.collapsed)
    left.distances[queryHits(left.hits)] <- 0
    right.distances[queryHits(right.hits)] <- 0
  }
  
  gr$leftDist <- left.distances
  gr$rightDist <- right.distances
  
  return(gr)
}

#' Function to calculate minimal distances to user defined ranges
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param userTrack A \code{\link{GRanges-class}} object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
getMinDist <- function(gr, userTrack) {
  distances <- range2rangeDistance(gr=gr, userTrack=userTrack, allow.overlap = TRUE)
  distances$leftDist[distances$leftDist == -1] <- distances$rightDist[distances$leftDist == -1]
  distances$rightDist[distances$rightDist == -1] <- distances$leftDist[distances$rightDist == -1]
  min.distances <- pmin(distances$leftDist, distances$rightDist)
  return(min.distances)
}