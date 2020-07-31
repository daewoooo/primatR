#' Synchronize Strand-seq read directionality
#'
#' This function aims to synchronize read directionality of a specified genomic region.
#' 
#' @param bamFiles A list of files that contains \code{\link{BreakPoint}} objects.
#' @param region A \code{\link{GRanges-class}} object to select reads from.
#' @inheritParams bam2GRanges
#' @return A \code{\link{GRanges-class}} object that reads synchronized by directionality.
#' @importFrom dplyr recode
#' @author David Porubsky
#' @export

synchronizeReadDirRegion <- function(bamFiles="", region = NULL, pairedEndReads = TRUE, min.mapq = 10, filterAltAlign = TRUE, alpha=0.1, genot.region.ends=TRUE) {
  #bams <- list.files(bamfolder, pattern = "\\.bam$", full.names = TRUE)
  
  if (is.null(region)) {
    "Please submit a genomic region. This function is not designed to operate in genome-wide settings!!!"
  }
  
  strand.sync <- list()
  data.merged <- GRanges()
  for (i in seq_along(bamFiles)) {
    bam <- bamFiles[i]
    filename <- basename(bam)
    message("Working on BAM: ", filename)
    
    ## Load BAM file
    data <- bam2GRanges(bamfile = bam, region = region, pairedEndReads = pairedEndReads, min.mapq =  min.mapq, filterAltAlign = filterAltAlign)

    ## Check n number of reads on either side of selected region have differing strand state by genotyping 25% of reads from each end of the given region
    if (genot.region.ends) {
      n.reads <- round(length(data) * 0.25)
      if (n.reads > 0) {
        left.end.reads <- data[1:n.reads]
        right.end.reads <- data[(1 + (length(data) - n.reads)):length(data)]
        left.end.probs <- countProb(minusCounts = length(left.end.reads[strand(left.end.reads) == '-']), 
                                    plusCounts = length(left.end.reads[strand(left.end.reads) == '+']), 
                                    alpha = alpha, log = TRUE)
        right.end.probs <- countProb(minusCounts = length(right.end.reads[strand(right.end.reads) == '-']), 
                                     plusCounts = length(right.end.reads[strand(right.end.reads) == '+']), 
                                     alpha = alpha, log = TRUE)
        ## Set data to NULL if the strand state on both sides of the given region is the same
        if (which.max(left.end.probs) == which.max(right.end.probs)) {
          data <- NULL
        }
      } else {
        data <- NULL
      }  
    }
    
    if (length(data) > 0) {
      ## Initialize master GRanges object to store all reads
      #if (i == init) {
      if (length(data.merged) == 0) {  
        read.dir <- as.vector(strand(data))
        read.dir <- dplyr::recode(read.dir, "+" = 1, "-" = 0)
        data$read.dir <- read.dir
        data.merged <- data
        strand.sync[[i]] <- data.frame(file=filename, dir="orig")
        init.idx <- i
      }
    
      ## Synchronize directionality and merge reads
      if (i > init.idx) {
        ## Recode read directionality into a binary code
        data$read.dir <- dplyr::recode(as.vector(strand(data)), "+" = 1, "-" = 0)
        ## Merge master GRanges with single cell GRanges and sort by position
        data.orig <- GenomicRanges::sort(c(data.merged, data), ignore.strand=TRUE)
        ## Calculate level of compression of a merge object
        data.orig.compression <- length(base::rle(data.orig$read.dir)$values)
        ## Switch directionality of single cell GRanges
        data$read.dir <- dplyr::recode(data$read.dir, "0" = 1, "1" = 0)
        strand(data) <- dplyr::recode(as.vector(strand(data)), "+" = "-", "-" = "+")
        ## Merge master GRanges with single cell GRanges (switched) and sort by position
        data.switch <- GenomicRanges::sort(c(data.merged, data), ignore.strand=TRUE)
        ## Calculate level of compression of a merge object
        data.switch.compression <- length(base::rle(data.switch$read.dir)$values)
        
        ## Compare level of compression between original ('orig') and switched ('switch') direction and decide if to keep directionality change
        if (data.orig.compression < data.switch.compression) {
          data.merged <- data.orig
          strand.sync[[i]] <- data.frame(file=filename, dir="orig")
        } else {
          data.merged <- data.switch
          strand.sync[[i]] <- data.frame(file=filename, dir="switch")
        }
      }
    }  
  }
  ## Make sure that max strand is Crick '+'
  max.strand <- names(which.max(BiocGenerics::table(strand(data.merged))))
  #switch strand if max directionality is Watson("-")
  if (max.strand == "-") {
    strand(data.merged) <- dplyr::recode(as.vector(strand(data.merged)), "+" = "-", "-" = "+")
    data.merged$read.dir <- dplyr::recode(data.merged$read.dir, "0" = 1, "1" = 0)
  }  
  ## Return final GRanges object with synchronized directionality
  return(data.merged)
}
