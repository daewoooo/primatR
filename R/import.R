#' Load structural variant calls
#'
#' @param file A file containing structural variant calls.
#' @param ID A metadata column to be appended.
#' @param remove.uncertain Set to \code{TRUE} if calls with question mark should be removed.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export

importCalls <- function(file=NULL, ID=NULL, remove.uncertain=TRUE) {
  calls.tab <- read.table(file = file, stringsAsFactors = FALSE, sep=" ", header = TRUE)
  ## Select only certain columns
  calls.tab <- calls.tab[,c('Chr.1','End.1','Start.2','SVclass','genoT')]
  ## Transform to GRanges
  calls.tab.gr <- GRanges(seqnames=calls.tab[,1], ranges=IRanges(start=calls.tab[,2], end=calls.tab[,3]), SVclass=calls.tab[,4], gen=calls.tab[,5])
  ## Remove uncertain calls with questionmark
  if (remove.uncertain) {
    calls.tab.gr <- calls.tab.gr[grep("\\?", calls.tab.gr$SVclass, invert=TRUE)]
  }
  ## Het is ‘+’, Hom is ‘-‘
  calls.tab.gr$gen[calls.tab.gr$gen == "+"] <- 'HET'
  calls.tab.gr$gen[calls.tab.gr$gen == "-"] <- 'HOM'
  ## ID column
  if (!is.null(ID)) {
    calls.tab.gr$ID <- as.character(ID)
  }
  return(calls.tab.gr)
}


#' Load reads from a composite file
#' 
#' This function loads Strand-seq reads stored as a RData object either all or from a selected genomic regions.
#'
#' @param compositeFile A file containing single-cell Strand-seq reads synchronized by directionality.
#' @param regions A \code{\link{GRanges-class}} object with genomic regions to select reads from.
#' @param ID A unique identifier to add as a metacolumn in returned \code{\link{GRanges-class}} object.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export

importReadsFromComposite <- function(compositeFile, regions=NULL, ID="") {
  data <- get(load(compositeFile))
  
  if (is.character(ID) & nchar(ID)>0) {
    data$ID <- ID
  } else {
    data$ID <- ""
  }
  
  if (!is.null(regions) & length(regions) > 0) {
    data <- IRanges::subsetByOverlaps(data, regions)
    hits <- GenomicRanges::findOverlaps(data, regions)
    data$region.ID <- rep(as.character(regions), table(subjectHits(hits)))
  }
  
  return(data)
}
