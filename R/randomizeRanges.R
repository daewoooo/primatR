#' Function to randomize user defined ranges along each chromosome
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param mask.gr A \code{\link{GRanges-class}} object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
randomizeRanges <- function(gr=NULL, mask.gr=NULL, bsgenome=NULL) {
  ## Add sequence lengths if missing
  if (any(is.na(seqlengths(gr)))) {
    ## Load BSgenome
    if (class(bsgenome) != 'BSgenome') {
      if (is.character(bsgenome)) {
        suppressPackageStartupMessages(library(bsgenome, character.only=T))
        bsgenome <- as.object(bsgenome) # replacing string by object
      }
    }
    seqlengths(gr) <- seqlengths(bsgenome)[names(seqlengths(gr))]
  }  
  
  shifted.grl <- GenomicRanges::GRangesList()
  for (i in 1:length(gr)) {
    range <- gr[i]
    chr.len <- seqlengths(range)[as.character(seqnames(range))]
    repeat {
      newStart <- round(runif(1, 1, chr.len))
      shift.gr <- GenomicRanges::GRanges(seqnames=seqnames(range), range=IRanges(start=newStart, width=width(range)))
      if (!is.null(mask.gr)) {
        if (GenomicRanges::countOverlaps(shift.gr, mask.gr) == 0) break
        if (end(shift.gr) < chr.len) break
      } else {
        if (end(shift.gr) < chr.len) break
      }
    }    
    shifted.grl[[i]] <- shift.gr
  }
  randShift.gr <- unlist(shifted.grl)
  return(randShift.gr)
}