#' Split Iso-seq alignments by aligned exons
#' 
#' This function takes aligned Iso-seq reads and splits every read into a sub-reads creads creating linked
#' read cloud of gene exons. 
#'
#' @param bamfile Bamfile with aligned Iso-seq reads.
#' @param out.format 
#' @param filename A filename to store exported sequences into.
#' @return A \code{NULL} object.
#' @importFrom Rsamtools scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments cigarRangesAlongReferenceSpace cigarToRleList
#' @author David Porubsky
#' @export
#' 
cutISSreads <- function(bamfile=NULL, out.format='fasta', filename=NULL) {
  ## Delete output file if it exists
  if (file.exists(filename)) {
    suppressMessages( file.remove(filename) )
  }
  ## Load BAM file
  message("Reading BAM file: ", basename(bamfile))
  data.raw <- GenomicAlignments::readGAlignments(bamfile, param=Rsamtools::ScanBamParam(what=c('cigar', 'mapq', 'flag', 'qname', 'seq'), flag=scanBamFlag(isDuplicate=FALSE)))
  frags <- as(data.raw, 'GRanges')
  ## Split aligned fragments by read name
  frags.grl <- split(frags, frags$qname)
  n.frags <- length(frags.grl)
  
  ## Go over every CCS iso-seq read
  counter <- 1000
  for (i in seq_along(frags.grl)) {
    if (i == counter) {
      message("    Processed [", i, "/", n.frags, "] fragments ...")
      counter <- counter + 1000
    }
    read.frags <- frags.grl[[i]]
    seq <- read.frags$seq[width(read.frags$seq) > 1] ## Get the raw sequences
    seq <- seq[which.max(width(seq))]

    ## Parse cigar string
    cigar.iranges <- GenomicAlignments::cigarRangesAlongQuerySpace(read.frags$cigar, flag = read.frags$flag)
    cigarRle <- GenomicAlignments::cigarToRleList(read.frags$cigar)
    ## Convert parsed cigar to GRanges
    cigar.gr <- GenomicRanges::GRanges(seqnames=as.character(seqnames(read.frags[1])), ranges=unlist(cigar.iranges))
    cigar.gr$cigar <- unlist(runValue(cigarRle))
    cigar.gr$qname <- factor(rep(read.frags$qname, lengths(cigar.iranges)), levels = unique(read.frags$qname))
    ## Get postions of splitted alignments
    cut.sites <- GenomicRanges::resize(cigar.gr[cigar.gr$cigar == 'N'], width = 1)
    seqlengths(cut.sites) <- nchar(seq)
    ## Get regions in-between the cut.sites
    seq.regions <- GenomicRanges::gaps(cut.sites)
    seq.regions <- seq.regions[strand(seq.regions) == '*']
    ## Get sequence corresponding to seq.regions
    seq.chunks <- BSgenome::Views(unlist(seq), ranges(seq.regions))
    seq.names <- paste0(unique(cigar.gr$qname), "_", seq(1:length(seq.chunks)))
    names(seq.chunks) <- seq.names
    seq.chunks <- as(seq.chunks, 'XStringSet')
    ## Write out cutted sequence
    Biostrings::writeXStringSet(x = seq.chunks, filepath = filename, append = TRUE, format = out.format, compress = TRUE)
  }
  message("DONE!!!")
}