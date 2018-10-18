#' Function to export raw inversion calls.
#' 
#' This function loads manually checked inversion calls and order them by position and exports corresponding UCSC formated file.
#'
#' @param infile A filepath to file that contains maually checked inversion calls.
#' @param outputfolder A path to a folder where the to save processed inversion calls.
#' @param index A unique index to distinguish inversion call from different individuals.
#' @return \code{NULL}
#' @author David Porubsky

exportINVcalls <- function(infile=NULL, outputfolder="./", index='') {
  inv.calls <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  
  ## Check if loaded table has all required columns [TODO]
  
  ## Create output directory if it doesn't exist
  if (!dir.exists(outputfolder)) {
    message("Creating directory ", outputfolder, " ...")
    dir.create(outputfolder)
  }
  
  ## Fill breakpoints for manually selected regions
  inv.calls[inv.calls$Start.1 == 0, 'Start.1'] <- inv.calls[inv.calls$Start.1 == 0, 'StartCI.1']
  inv.calls[inv.calls$End.1 == 0, 'End.1'] <- inv.calls[inv.calls$End.1 == 0, 'EndCI.1']
  inv.calls[inv.calls$Start.2 == 0, 'Start.2'] <- inv.calls[inv.calls$Start.2 == 0, 'StartCI.2']
  inv.calls[inv.calls$End.2 == 0, 'End.2'] <- inv.calls[inv.calls$End.2 == 0, 'EndCI.2']
  
  inv.calls.gr <- GenomicRanges::GRanges(
    seqnames = inv.calls$Chr.1, 
    ranges=IRanges(start=inv.calls$Start.1, end=inv.calls$End.2), 
    SVclass=inv.calls$SVclass, 
    genoT=inv.calls$genoT
  )
  
  ## Order inversion calls by genomic position
  region.ord <- order(inv.calls.gr)
  inv.calls.gr <- inv.calls.gr[region.ord]
  
  inv.calls.ordered <- inv.calls[region.ord,]
  inv.calls.ordered <- inv.calls.ordered[,-c(5,6,13,14)]
  
  if (nchar(index) == 0) { index <- 'unknown' }
  
  destination <- file.path(outputfolder, paste0(index, "_INVcalls.ordered.txt"))
  write.table(inv.calls.ordered, file = destination, quote = FALSE, row.names = FALSE)
  
  ## Filter only simple inversions [OPTIONAL]
  #inv.calls.ordered.INVonly <- inv.calls.ordered[inv.calls.ordered$SVclass == 'INV']
  #inv.calls.ordered.INVonly <- inv.calls.ordered[inv.calls.ordered$SVclass == 'INV',]
  
  ## Export UCSC bed formated files
  inv.calls.ucsc <- inv.calls.ordered
  inv.calls.ucsc <- inv.calls.ucsc[,c('Chr.1', 'Start.2', 'End.1', 'genoT', 'SVclass')]
  inv.calls.ucsc$score <- 0
  ## Color genotypes by strand 'HET' == '+' && 'HOM' == '-'
  inv.calls.ucsc <- inv.calls.ucsc[,c('Chr.1', 'Start.2', 'End.1', 'SVclass', 'score', 'genoT')]
  
  ## Write to a file
  ucsc.header <- paste0('track name=', index,' ROIs description=', index, '_Bed_of_inverted_regions visibility=dense colorByStrand=\"28,144,153 221,28,119\"')
  destination <- file.path(outputfolder, paste0(index, "_INVcalls.bed"))
  write.table(ucsc.header, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=FALSE, sep='\t')
  write.table(inv.calls.ucsc, file=destination, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep='\t')
}