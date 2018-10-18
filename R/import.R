#' Load inversion calls
#'
#' @param file A file containing inversion calls.
#' @param ID A metadata column to be appended.
#' @param remove.uncertain Set to \code{TRUE} if calls with question mark should be removed.
#' @return A \code{\link{GRanges-class}} object containing all inversion calls.
#' @author David Porubsky

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
