#' Function to randomize user defined ranges along each chromosome
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param mask.gr A \code{\link{GRanges-class}} object.
#' @param bsgenome A \code{\link{GBSgenome-class}} object to get sequence lengths of submitted 'gr' object.
#' @param fai A FASTA index to get sequence lengths of submitted 'gr' object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
randomizeRanges <- function(gr=NULL, mask.gr=NULL, bsgenome=NULL, fai=NULL) {
  ## Add sequence lengths if missing
  if (any(is.na(seqlengths(gr)))) {
    if (!is.null(bsgenome)) {
      ## Load BSgenome
      if (class(bsgenome) != 'BSgenome') {
        if (is.character(bsgenome)) {
          suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
          bsgenome <- as.object(bsgenome) # replacing string by object
        }
      }
      seqlengths(gr) <- GenomeInfoDb::seqlengths(bsgenome)[names(seqlengths(gr))]
    } else if (!is.null(fai)) {
      ## Get seqeunce lengths from fasta index
      fai.tab <- utils::read.table(fai)
      fai.tab <- fai.tab[order(fai.tab$V2, decreasing = TRUE),]
      chrom.lengths <- fai.tab$V2
      names(chrom.lengths) <- fai.tab$V1
      seqlengths(gr) <- chrom.lengths[names(seqlengths(gr))]
    } else {
      stop("Submitted 'GRanges' object doesn't contain sequence lengths. Sequence lengths can be obtained from 'bsgenome' object or from fasta index '.fai'!!!")
    } 
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
