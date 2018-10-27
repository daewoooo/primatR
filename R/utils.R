#' Insert chromosome for in case it's missing
#' 
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} object with an additional metadata column containing chromosome name with 'chr'. 
insertchr <- function(gr) {
  mask <- which(!grepl('chr', seqnames(gr)))
  mcols(gr)$chromosome <- as.character(seqnames(gr))
  mcols(gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(gr)$chromosome[mask])
  mcols(gr)$chromosome <- as.factor(mcols(gr)$chromosome)
  return(gr)
}