#' Detect trans-chromosomal mappings of Iso-seq reads
#'
#' This function reads BAM file into a \code{\link{GRanges}} object and searches for reads that
#' map do different chromosomal locations. Links supported by user defined number of reads are 
#' reported and points to a putative gene fusions.
#'
#' @param bamfile Bamfile with aligned reads.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param standardChroms If set to \code{TRUE} only standard chromosomes (1-22, X & Y) will be kept.
#' @param min.links Minimum number of supporting Iso-seq reads in order to report a putative gene fusion.
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @return A \code{\link{GRanges}} object reporting all found links.
#' @author David Porubsky
#' @export
#' 
getTransChromFusions <- function(bamfile=NULL, min.mapq=10, standardChroms=TRUE, min.links=5) {
  message("Reading BAM file: ", basename(bamfile))
  ## Load BAM file
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignments(bamfile, param=Rsamtools::ScanBamParam(what=c('mapq','qname'), flag=scanBamFlag(isDuplicate=FALSE))) )
  frags <- as(data.raw, 'GRanges')
  ## Keep only standard chromosomes
  if (standardChroms) {
    frags <- GenomeInfoDb::keepStandardChromosomes(frags, pruning.mode = 'coarse')
  }  
  ## Filter by mapping quality
  if (!is.null(min.mapq)) {
    if (any(is.na(mcols(frags)$mapq))) {
      warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
      mcols(frags)$mapq[is.na(mcols(frags)$mapq)] <- -1
    }
    frags <- frags[mcols(frags)$mapq >= min.mapq]
  }
  ## Split aligned fragments by read name
  frags.grl <- split(frags, frags$qname)
  pos.per.frag <- lengths(frags.grl)
  ## Select reads that maps to multiple locations
  frag.multi.pos <- which(pos.per.frag > 1)
  frags.grl <- frags.grl[frag.multi.pos]
  ## Keep reads that maps to different chromosomes
  filt <- sapply(frags.grl, function(x) length(unique(as.character(seqnames(x)))) > 1)
  frags.grl <- frags.grl[filt]
  ## Get region links
  frags.gr <- unlist(frags.grl, use.names = FALSE)
  ## Get non-overlappign genomic regions
  regions.gr <- frags.gr
  strand(regions.gr) <- "*"
  regions.gr <- GenomicRanges::reduce(regions.gr)
  regions.gr$idx <- 1:length(regions.gr)
  ## Create empty data matrix to store links
  links.m <- matrix(0L, nrow = length(regions.gr), ncol = length(regions.gr))
  ## Process every read and record link indices
  for (i in seq_along(frags.grl)) {
    frag <- frags.grl[[i]]
    ## Skip read that map to more than 2 genomic locations
    if (length(runValue(seqnames(frag))) > 2) {
      next
    }
    ## Skip read that map to more than 2 genomic locations within a chromosome
    if (any(runLength(seqnames(frag)) > 1)) {
      next
    }
    hits <- findOverlaps(frag, regions.gr)
    region.idx <- sort(unique(subjectHits(hits)))
    #message("iter:",i , " ", length(region.idx))
    ## Increment score for a given pair of regions (indices)
    links.m[region.idx[1], region.idx[2]] <- links.m[region.idx[1], region.idx[2]] + 1
  }
  ## Get coordinates of candidate fusions with 5 and more supporting Iso-seq reads
  fusion.idx <- which(links.m >= min.links, arr.ind = TRUE)
  scores <- diag(links.m[fusion.idx[,1], fusion.idx[,2]])
  
  fusion.links <- regions.gr[fusion.idx[,1], 0]
  fusion.links$to.gr <- regions.gr[fusion.idx[,2], 0]
  fusion.links$scores <- scores
  ## Export final results
  return(fusion.links)
}  


