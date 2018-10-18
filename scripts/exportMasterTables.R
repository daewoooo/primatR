## Export master table for each great ape ##
############################################
## Chimp ##
strandS.dorien.gr.export <- strandS.dorien.gr[strandS.dorien.gr$SVclass == 'INV']
strandS.dorien.gr.export <- gr2grOverlap(query.gr = strandS.dorien.gr.export, subject.gr = seg.dup.gr)

strandS.dorien.gr.export <-  getReciprocalOverlaps(gr1 = strandS.dorien.gr.export, gr2 = catacchio.calls.panTro.gr, thresh = 0, return = 'query')
strandS.dorien.gr.export$Catacchio.perc.overlap <- strandS.dorien.gr.export$perc.overlap
strandS.dorien.gr.export$Catacchio.range <- strandS.dorien.gr.export$toGR
strandS.dorien.gr.export$Catacchio.valid <- strandS.dorien.gr.export$perc.overlap >= 50

strandS.dorien.gr.export <-  getReciprocalOverlaps(gr1 = strandS.dorien.gr.export, gr2 = zev.callsValid.chimp.gr, thresh = 0, return = 'query')
strandS.dorien.gr.export$Zev.perc.overlap <- strandS.dorien.gr.export$perc.overlap
strandS.dorien.gr.export$Zev.range <- strandS.dorien.gr.export$toGR
strandS.dorien.gr.export$Zev.valid <- strandS.dorien.gr.export$perc.overlap >= 50

strandS.dorien.gr.export <-  getReciprocalOverlaps(gr1 = strandS.dorien.gr.export, gr2 = feuk.callsChimp.gr, thresh = 0, return = 'query')
strandS.dorien.gr.export$Feuk.perc.overlap <- strandS.dorien.gr.export$perc.overlap
strandS.dorien.gr.export$Feuk.range <- strandS.dorien.gr.export$toGR
strandS.dorien.gr.export$Feuk.valid <- strandS.dorien.gr.export$perc.overlap >= 50

strandS.dorien.gr.export <-  getReciprocalOverlaps(gr1 = strandS.dorien.gr.export, gr2 = BioNano.chimp.gr, thresh = 0, return = 'query')
strandS.dorien.gr.export$BioN.perc.overlap <- strandS.dorien.gr.export$perc.overlap
strandS.dorien.gr.export$BioN.range <- strandS.dorien.gr.export$toGR
strandS.dorien.gr.export$BioN.valid <- strandS.dorien.gr.export$perc.overlap >= 50

remove.cols <- which( names(mcols(strandS.dorien.gr.export)) %in% c('perc.overlap','toGR','idx') )
strandS.dorien.gr.export <- strandS.dorien.gr.export[,-remove.cols]
strandS.dorien.df.export <- as(strandS.dorien.gr.export, 'data.frame')
write.table(strandS.dorien.df.export , file = "/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/chimp_master_table.csv", quote = FALSE, sep = ",", row.names = FALSE)


## Bonobo ##
strandS.ulindi.gr.export <- strandS.ulindi.gr[strandS.ulindi.gr$SVclass == 'INV']
strandS.ulindi.gr.export <- gr2grOverlap(query.gr = strandS.ulindi.gr.export, subject.gr = seg.dup.gr)

strandS.ulindi.gr.export <-  getReciprocalOverlaps(gr1 = strandS.ulindi.gr.export, gr2 = BioNano.bonobo.gr, thresh = 0, return = 'query')
strandS.ulindi.gr.export$BioN.perc.overlap <- strandS.ulindi.gr.export$perc.overlap
strandS.ulindi.gr.export$BioN.range <- strandS.ulindi.gr.export$toGR
strandS.ulindi.gr.export$BioN.valid <- strandS.ulindi.gr.export$perc.overlap >= 50

remove.cols <- which( names(mcols(strandS.ulindi.gr.export)) %in% c('perc.overlap','toGR','idx') )
strandS.ulindi.gr.export <- strandS.ulindi.gr.export[,-remove.cols]
strandS.ulindi.df.export <- as(strandS.ulindi.gr.export, 'data.frame')
write.table(strandS.ulindi.df.export , file = "/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/Final_callset/bonobo_master_table.csv", quote = FALSE, sep = ",", row.names = FALSE)


## Orangutan ##
strandS.ppy.gr.export <- strandS.ppy.gr[strandS.ppy.gr$SVclass == 'INV']
strandS.ppy.gr.export <- gr2grOverlap(query.gr = strandS.ppy.gr.export , subject.gr = seg.dup.gr)

strandS.ppy.gr.export <-  getReciprocalOverlaps(gr1 = strandS.ppy.gr.export, gr2 = catacchio.calls.ppy.gr, thresh = 0, return = 'query')
strandS.ppy.gr.export$Catacchio.perc.overlap <- strandS.ppy.gr.export$perc.overlap
strandS.ppy.gr.export$Catacchio.range <- strandS.ppy.gr.export$toGR
strandS.ppy.gr.export$Catacchio.valid <- strandS.ppy.gr.export$perc.overlap >= 50

strandS.ppy.gr.export <-  getReciprocalOverlaps(gr1 = strandS.ppy.gr.export, gr2 = zev.callsValid.ppy.gr, thresh = 0, return = 'query')
strandS.ppy.gr.export$Zev.perc.overlap <- strandS.ppy.gr.export$perc.overlap
strandS.ppy.gr.export$Zev.range <- strandS.ppy.gr.export$toGR
strandS.ppy.gr.export$Zev.valid <- strandS.ppy.gr.export$perc.overlap >= 50

strandS.ppy.gr.export <-  getReciprocalOverlaps(gr1 = strandS.ppy.gr.export, gr2 = BioNano.ppy.gr, thresh = 0, return = 'query')
strandS.ppy.gr.export$BioN.perc.overlap <- strandS.ppy.gr.export$perc.overlap
strandS.ppy.gr.export$BioN.range <- strandS.ppy.gr.export$toGR
strandS.ppy.gr.export$BioN.valid <- strandS.ppy.gr.export$perc.overlap >= 50

remove.cols <- which( names(mcols(strandS.ppy.gr.export)) %in% c('perc.overlap','toGR','idx') )
strandS.ppy.gr.export <- strandS.ppy.gr.export[,-remove.cols]
strandS.ppy.gr.export <- as(strandS.ppy.gr.export, 'data.frame')
write.table(strandS.ppy.gr.export , file = "/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_500KbBin_withAlt/orangutan_master_table.csv", quote = FALSE, sep = ",", row.names = FALSE)