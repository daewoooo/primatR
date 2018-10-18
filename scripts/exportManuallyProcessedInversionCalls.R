## Export manually processed inversion calls ##
###############################################

library(primatR)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Chimpanzee/Dorien_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt_old/syncReads_Multireads_chimp_breakPoints_checked.bed", 
               outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls",
               index = "chimpanzee"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Bonobo/Ulindi_Bams_GRCh38/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_bonobo_breakPoints_checked.bed", 
               outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls",
               index = "bonobo"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Gorilla/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_gorilla_breakPoints_checked.bed", 
               outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls",
               index = "gorilla"
)

exportINVcalls(infile = "/home/porubsky/WORK/Great_apes/Orangutan/Ashley_bam/selected/BreakpointR_results_1MbBin_filtAlt/syncReads_Multireads_orangutan_breakPoints_checked.bed", 
               outputfolder = "/home/porubsky/WORK/Great_apes/Final_INV_calls",
               index = "orangutan"
)