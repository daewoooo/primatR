#' Synchronize Strand-seq read directionality
#'
#' This function aims to synchronize read directionality of a specified genomic region.
#' 
#' @param bamfiles A list of files that contains \code{\link{BreakPoint}} objects.
#' @param region A \code{\link{GRanges-class}} object to select reads from.
#' @param genot.frac.reads A fraction of the total number of reads to be genotype from both or a single end 
#' of a given 'region' based on 'genot.region.ends'. 
#' @param genot.region.ends Define if 'left' or 'right' of 'both ends of a given region shall be genotyped.
#' @param desired.genot Select desired genotypes: 'ww', 'cc', 'wc' to be allowed on selected region ends. If 
#' 'genot.region.ends' = 'both', both ends have to have the same genotype.
#' @inheritParams bam2GRanges
#' @inheritParams countProb
#' @return A \code{\link{GRanges-class}} object that reads synchronized by directionality.
#' @importFrom dplyr recode
#' @importFrom BiocGenerics table
#' @author David Porubsky
#' @export

synchronizeReadDirRegion <- function(bamfiles="", region = NULL, pairedEndReads = TRUE, min.mapq = 10, filterAltAlign = TRUE, alpha=0.1, genot.frac.reads=0.25, genot.region.ends='right', desired.genot=c('cc','ww')) {
  message("    Synchronizing read directionality for region: ", as.character(region), " ...", appendLF=F); ptm <- proc.time()
  
  if (is.null(region)) {
    stop("Please submit a genomic region as GRanges object. This function is not designed to operate in genome-wide settings!!!")
  }
  
  strand.sync <- list()
  data.merged <- GRanges()
  for (i in seq_along(bamfiles)) {
    bam <- bamfiles[i]
    filename <- basename(bam)
    #message("Working on BAM: ", filename)
    
    ## Load BAM file
    data <- bam2GRanges(bamfile = bam, region = region, pairedEndReads = pairedEndReads, min.mapq =  min.mapq, filterAltAlign = filterAltAlign)
    
    ## Set desired genotypes at region ends
    if (is.character(desired.genot)) {
      avail.genot <- list('ww' = 1, 'cc' = 2, 'wc' = 3)
      search.genot <- avail.genot[grep(names(avail.genot), pattern = paste(desired.genot, collapse = "|"), ignore.case = TRUE)]
      search.genot <- unlist(search.genot)
    } else {
      search.genot <- c(1:3)
    }
    
    ## Check n number of reads on either side of selected region have differing strand state by genotyping 'genot.frac.reads' of reads from each end of the given region
    n.reads <- round(length(data) * genot.frac.reads)
    if (genot.region.ends == 'right') {
      if (n.reads > 0) {
        right.end.reads <- data[(1 + (length(data) - n.reads)):length(data)]
        right.end.probs <- countProb(minusCounts = length(right.end.reads[strand(right.end.reads) == '-']), 
                                     plusCounts = length(right.end.reads[strand(right.end.reads) == '+']), 
                                     alpha = alpha, log = TRUE)
        if (!which.max(right.end.probs) %in% search.genot) {
          data <- NULL
        }
      } else {
        data <- NULL
      }  
    } else if (genot.region.ends == 'left') {
      left.end.reads <- data[1:n.reads]
      left.end.probs <- countProb(minusCounts = length(left.end.reads[strand(left.end.reads) == '-']), 
                                  plusCounts = length(left.end.reads[strand(left.end.reads) == '+']), 
                                  alpha = alpha, log = TRUE)
      if (!which.max(right.end.probs) %in% search.genot) {
        data <- NULL
      }
    } else if (genot.region.ends == 'both') {
      if (n.reads > 0) {
        right.end.reads <- data[(1 + (length(data) - n.reads)):length(data)]
        left.end.reads <- data[1:n.reads]
        right.end.probs <- countProb(minusCounts = length(right.end.reads[strand(right.end.reads) == '-']), 
                                     plusCounts = length(right.end.reads[strand(right.end.reads) == '+']), 
                                     alpha = alpha, log = TRUE)
        left.end.probs <- countProb(minusCounts = length(left.end.reads[strand(left.end.reads) == '-']), 
                                    plusCounts = length(left.end.reads[strand(left.end.reads) == '+']), 
                                    alpha = alpha, log = TRUE)
        if (which.max(right.end.probs) != which.max(left.end.probs)) {
          data <- NULL
        }
        if (!which.max(right.end.probs) %in% search.genot | !which.max(left.end.probs) %in% search.genot) {
          data <- NULL
        }
      } else {
        data <- NULL
      }   
    } else {
      data <- NULL
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
  ## Report time
  time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
  ## Return final GRanges object with synchronized directionality
  return(data.merged)
}