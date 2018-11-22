## Load required libraries
library(GenomicAlignments)
library(RColorBrewer)
#library(rbamtools)
library(ggplot2)
library(scales)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyr)

## Set directories
data.dir <- "/home/porubsky/WORK/Great_apes/InvertedDups_analysis/TEST2/"
bamfiles <- list.files("/home/porubsky/WORK/Great_apes/InvertedDups_analysis/TEST2/", pattern = "\\.bam$", full.names = TRUE)

## Set parameters
mapq <- 60
filt.flag <- 3328
chromosomes <- paste0('chr', c(1:22, 'X'))
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

## Load segDup Track
seg.dup <- read.table("/home/porubsky/WORK/Great_apes/CompareWithPreviousStudies/segDupTrackUCSC_hg38.bed")
seg.dup.gr <- GRanges(seqnames=seg.dup$V1, ranges=IRanges(start=seg.dup$V2, end=seg.dup$V3), fracMatch=seg.dup$V4)

fragments <- GRangesList()
for (i in seq_along(bamfiles)) {
  bam <- bamfiles[[i]]
  filename <- basename(bam)
  message("Processing BAM: ", filename)
  
  ## Load read pairs
  suppressWarnings( data.raw <- GenomicAlignments::readGAlignmentPairs(bam, param=Rsamtools::ScanBamParam(what=c('mapq', 'flag', 'qname'))) )
  ## Convert to GRanges
  data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
  data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
  
  ## Remove positions overlapping with SDs
  hits.first <- findOverlaps(data.first, seg.dup.gr)
  hits.last <- findOverlaps(data.last, seg.dup.gr)
  mask <- unique(sort(c(queryHits(hits.first), queryHits(hits.last))))
  data.first <- data.first[-mask]
  data.last <- data.last[-mask]
  
  ## Filter by mapq
  if (mapq > 0) {
    mask <- data.first$mapq >= mapq & data.last$mapq >= mapq
    data.first <- data.first[mask]
    data.last <- data.last[mask]
  }

  ## Filter by flag
  if (filt.flag > 0) {
    bit.flag.first <- bitwAnd(filt.flag, data.first$flag)
    bit.flag.last <- bitwAnd(filt.flag, data.last$flag)
    mask <- bit.flag.first == 0 & bit.flag.last == 0
    data.first <- data.first[mask]
    data.last <- data.last[mask]
  } 
  
  ## Export fragments
  data <- data.first
  data$to.gr <- data.last
  data$ID <- sapply(data$qname, function(x) strsplit(x, "__")[[1]][2])
  fragments[[i]] <- data
}
  
## Merge all data for plotting
data.plt <- unlist(fragments, use.names = FALSE)
assessed.breaks <- length(unique(data.plt$ID))

## TODO filter random links here???
binned.gr <- makeBins(bsgenome = bsgenome, chromosomes = chromosomes, binsize = 100000, stepsize = 50000)
from.gr <- data.plt[,'qname']
to.gr <- data.plt$to.gr[,'qname']
all.gr <- c(from.gr, to.gr)

hits <- findOverlaps(binned.gr, all.gr)
readIDs.per.overlap <- split(all.gr$qname[subjectHits(hits)], queryHits(hits))
unique.readIDs.per.overlap <- sapply(readIDs.per.overlap, countUniqueReadIDs)
## Get bins with minimal number of unique reads
bin.idx <- as.numeric( names(unique.readIDs.per.overlap)[unique.readIDs.per.overlap >= 10] )
select.regions <- binned.gr[bin.idx]
## Collapse overlapping ranges
filt.regions <- reduce(select.regions)

#binned.gr$from.gr.counts <- countOverlaps(binned.gr, from.gr)
#binned.gr$to.gr.counts <- countOverlaps(binned.gr, to.gr)
#binned.gr$total.counts <- binned.gr$from.gr.counts + binned.gr$to.gr.counts
#filt.regions <- binned.gr[binned.gr$total.counts > 50]
## Find overlapping ranges
from.gr.hits <- findOverlaps(from.gr, filt.regions)
to.gr.hits <- findOverlaps(to.gr, filt.regions)
filt.idx <- intersect(queryHits(from.gr.hits), queryHits(to.gr.hits))
data.plt <- data.plt[filt.idx]

## Prepare data.frame for plotting
plt.df <- as.data.frame(data.plt)
plt.df$break.chr <- sapply(plt.df$ID, function(x) gsub(strsplit(x, "-")[[1]][1], pattern = 'chr', replacement = ''))

## Add x_coord by selecting larger number being the start of the link
plt.df$x <- pmax(plt.df$start, plt.df$to.gr.start)
plt.df$xend <- pmin(plt.df$start, plt.df$to.gr.start)
## Remove links with the same start and end position
plt.df$link <- paste0(plt.df$x, "_", plt.df$xend)
plt.df <- plt.df[!duplicated(plt.df$link),]
## Remove links with the same start and end position
plt.df <- plt.df[plt.df$x != plt.df$xend,]

## Filter only standard chromosomes [1:22, X]
mask <- plt.df$seqnames %in% chromosomes & plt.df$to.gr.seqnames %in% chromosomes
plt.df <- plt.df[mask,]

## Impose chromosome levels to observed data
plt.df$y <- as.character(plt.df$seqnames)
plt.df$yend <- as.character(plt.df$to.gr.seqnames)
plt.df$y <- match(plt.df$y, as.character(chromosomes))
plt.df$yend <- match(plt.df$yend, as.character(chromosomes))

## Make sure that at least one break resides on chromosome of origin
#mask <- plt.df$break.chr != plt.df$y & plt.df$break.chr != plt.df$yend
plt.df <- plt.df %>% filter(break.chr == y | break.chr == yend)

## Add links between chromosomes
chr.from <- pmax(plt.df$y, plt.df$yend)
chr.to <- pmin(plt.df$y, plt.df$yend)
plt.df$chr.link <- paste0(chr.from, "_", chr.to)

## Get number of unique reads per link
qnames.list <- split(plt.df$qname, plt.df$chr.link)
qnames.counts.perLink <- sapply(qnames.list , countUniqueReadIDs)
signif.links <- names(qnames.counts.perLink)[qnames.counts.perLink >= 10]
plt.df.filt <- plt.df[plt.df$chr.link %in% signif.links,]

## Get number of intra- versus inter-chromosomal links
total.breaks <- length(unique(plt.df.filt$ID))
breaks.df <- split(plt.df.filt, plt.df.filt$ID)
intra.links <- 0
inter.links <- 0
for (i in seq_along(breaks.df)) {
  break.df <- breaks.df[[i]]
  links <- unique(break.df$chr.link)
  for (j in seq_along(links)) {
    link <- links[j]
    chroms <- strsplit(link, "_")[[1]]
  
    if (chroms[1] == chroms[2]) {
      intra.links <- intra.links + 1
    } else {
      inter.links <- inter.links + 1
    }
  }  
}

## Prepare genome-wide ideogram
seq.len <- seqlengths(bsgenome)[chromosomes]
ideo.df <- data.frame(seqnames=names(seq.len), length=seq.len)
ideo.df$seqnames <- factor(ideo.df$seqnames, levels=chromosomes)
ideo.df$levels <- 1:length(seq.len)
ideo <- ggplot(ideo.df) + geom_linerange(aes(x=levels, ymin=0, ymax=length), size=1, color='black')

## Add links to the ideogram
plt <- ideo +
  geom_curve(data=plt.df.filt, aes(x = y, y = start, xend = yend, yend = to.gr.start), color="red", curvature = 0.25) + 
  scale_y_continuous(labels = comma) +
  scale_x_continuous(breaks = ideo.df$levels, labels = ideo.df$seqnames) +
  guides(colour = guide_legend(override.aes = list(alpha=1))) +
  theme_bw() +
  ylab("Genomic position (bp)") +
  xlab("Chromosome")

## Add some statistics
plt <- plt + annotate(geom="text", x=23, y=250000000, label=paste0("Assessed breaks: ", assesse.breaks, "\n"), color="black", hjust=1)
plt <- plt + annotate(geom="text", x=23, y=250000000-10000000, label=paste0("Intra-chromosomal breaks: ", intra.links, "\n"), color="black", hjust=1)
plt <- plt + annotate(geom="text", x=23, y=250000000-20000000, label=paste0("Inter-chromosomal breaks: ", inter.links, "\n"), color="black", hjust=1)


##########################################################################################################
## Functions

makeBins <- function(bsgenome, chromosomes, binsize=100000, stepsize=binsize/2) {
  bins <- GRangesList()
  chr.lengths <- seqlengths(bsgenome)[chromosomes]
  seqlevels(bins) <- chromosomes
  for (i in seq_along(chr.lengths)) {
    chr.len <- chr.lengths[i]
    
    bin.starts <- seq(from = 1, to = chr.len-binsize, by = stepsize)
    bin.ends <- seq(from = binsize, to = chr.len, by = stepsize)

    chr.bins <- GRanges(seqnames=names(chr.len), ranges=IRanges(start=bin.starts, end=bin.ends))
    bins[[i]] <- chr.bins
  }
  bins <- unlist(bins, use.names = FALSE)
  seqlengths(bins) <- chr.lengths
  return(bins)
}

countUniqueReadIDs <- function(readIDs) {
  read.id <- sapply(readIDs, function(x) strsplit(x, "__")[[1]][1])
  unique.ids <- unique(read.id)
  return(length(unique.ids))
} 
