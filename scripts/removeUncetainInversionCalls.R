library(primatR)

## Load Strand-seq calls ##
###########################
## Strand-seq calls Chimp
strandS.dorien <- importCalls(file = "/home/porubsky/WORK/Great_apes/Final_INV_calls/chimpanzee_INVcalls.ordered.txt", ID = "Chimp", remove.uncertain = FALSE)
## Strand-seq calls Bonobo
strandS.ulindi <- importCalls(file = "/home/porubsky/WORK/Great_apes/Final_INV_calls/bonobo_INVcalls.ordered.txt", ID = "Bonobo", remove.uncertain = FALSE)
## Strand-seq calls Gorilla
strandS.ggo9 <- importCalls(file = "/home/porubsky/WORK/Great_apes/Final_INV_calls/gorilla_INVcalls.ordered.txt", ID = "GGO9", remove.uncertain = FALSE)
## Strand-seq calls Orangutan
strandS.ppy <- importCalls(file = "/home/porubsky/WORK/Great_apes/Final_INV_calls/orangutan_INVcalls.ordered.txt", ID = "PPY", remove.uncertain = FALSE)

## Retain uncertain calls based on PacBio calls (Zev et al. 2018) ##
####################################################################
## Published calls for chimpanzee, gorilla and orangutan
smartieSV.calls <- read.table(file = "/home/porubsky/WORK/Great_apes/CompareWithPreviousStudies/Zev_2018/aar6343_TableS12_2.csv", sep = ",", header = TRUE)
smartieSV.calls.ID <- smartieSV.calls$ID
smartieSV.calls.ID <- strsplit(as.character(smartieSV.calls.ID), "_")
smartieSV.calls.ID <- sapply(smartieSV.calls.ID, function(x) x[4])
smartieSV.calls.gr <- GRanges(seqnames=smartieSV.calls$GRCh38_seqid, ranges=IRanges(start=smartieSV.calls$start, end=smartieSV.calls$end), ID=smartieSV.calls.ID)
smartieSV.calls.grl <- split(smartieSV.calls.gr, smartieSV.calls.gr$ID)

strandS.dorien.uncertain <- strandS.dorien[grep("\\?", strandS.dorien$SVclass)]
dorien.uncertainCalls <- getReciprocalOverlaps(strandS.dorien.uncertain, smartieSV.calls.grl[['chimpanzee']], thresh = 0, return = 'query')

strandS.ulindi.uncertain <- strandS.ulindi[grep("\\?", strandS.ulindi$SVclass)] ## no SV calls

strandS.ggo9.uncertain <- strandS.ggo9[grep("\\?", strandS.ggo9$SVclass)]
ggo9.uncertainCalls <- getReciprocalOverlaps(strandS.ggo9.uncertain, smartieSV.calls.grl[['gorilla']], thresh = 0, return = 'query')

strandS.ppy.uncertain <- strandS.ppy[grep("\\?", strandS.ppy$SVclass)]
ppy.uncertainCalls <- getReciprocalOverlaps(strandS.ppy.uncertain, smartieSV.calls.grl[['orangutan']], thresh = 0, return = 'query')

validated.uncertainCalls <- c(dorien.uncertainCalls, ggo9.uncertainCalls, ppy.uncertainCalls)
validated.uncertainCalls <- validated.uncertainCalls[validated.uncertainCalls$perc.overlap>0] #Remove uncertain calls with no support in PacBio data
validated.uncertainCalls$SVclass <- gsub(validated.uncertainCalls$SVclass, pattern = "\\?", replacement = "") #Remove remaining unvalidated uncertain calls


all.gr <- all.gr[grep("\\?", all.gr$SVclass, invert = TRUE)]