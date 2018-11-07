#' Synchronize Strand-seq read directionality
#'
#' This function aims to synchronize read directionality of a specified genomic region.
#' 
#' @param bamfolder A list of files that contains \code{\link{BreakPoint}} objects.
#' @param region A segment size to be collapsed with neighbouring segments.
#' @inheritParams bam2GRanges
#' @return A \code{\link{GRanges-class}} object that reads synchronized by directionality.
#' @importFrom dplyr recode
#' @author David Porubsky
#' @export

synchronizeReadDirRegion <- function(bamfolder = "", region = NULL, pairedEndReads = TRUE, min.mapq = 10, filterAltAlign = TRUE) {
  bams <- list.files(bamfolder, pattern = "\\.bam$", full.names = TRUE)
  
  if (is.null(region)) {
    "Please submit genomic region. This function has not been tested in genome-wide settings!!!"
  }
  
  for (i in seq_along(bams)) {
    bam <- bams[i]
    filename <- basename(bam)
    message("Working on BAM: ", filename)
    
    ## Load BAM file
    data <- bam2GRanges(bamfile = bam, region = region, pairedEndReads = pairedEndReads, min.mapq =  min.mapq, filterAltAlign = filterAltAlign)

    ## Initialize master GRanges object to store all reads
    if (i == 1) {
      read.dir <- as.vector(strand(data))
      read.dir <- dplyr::recode(read.dir, "+" = 1, "-" = 0)
      data$read.dir <- read.dir
      data.merged <- data
      strand.sync[[i]] <- data.frame(file=filename, dir="orig")
    }
    
    ## Synchronize directionality and merge reads
    if (i > 1) {
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
