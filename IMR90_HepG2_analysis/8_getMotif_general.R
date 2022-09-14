###########################################################################################
########################## Get motifs for ChIPanalyser - general ##########################
###########################################################################################


# This whole script has been generalised apart from the lines in bash e.g. the ones below and at the end.
# screen, qrsh, cd /, go to IMR90/preProcessing folder
ls *treat_pileup.bdg > IMR90SEbdg.txt
mv /storage/st20d/HiC_Kc167/alessandra/IMR90/preProcessing/IMR90SEbdg.txt /storage/st20d/HiC_Kc167/alessandra/IMR90/ChIPanalyser/motifs/
cd ChIPanalyser/motifs/
R3.5.0

library(MotifDb)

cell_line <- "IMR90"

filename <- read.csv(paste0(cell_line,"SEbdg.txt"), sep ="\t", head = F)
TFnames <- gsub(paste0("-human_",cell_line,"_SE_trim-peaks_treat_pileup.bdg"), "", filename$V1) # get TF names


## Get motifs for each TF
for(i in TFnames){
  
  if(length(query(query(query(MotifDb, "hsapiens",
                              notStrings = c("Mmusculus", "rrattus")), i), "jaspar2018"))!= 0){

    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "jaspar2018")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    ##get the source of the PFM
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "jaspar2018"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10"))!= 0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10")
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02")) !=0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "jolma2013")) !=0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "jolma2013")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(MotifDb, "hsapiens",
                              notStrings=c("Mmusculus", "rrattus")), i))!= 0){
    PFM <- query(query(MotifDb, "hsapiens", notStrings=c("Mmusculus", "rrattus")), i)[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(MotifDb, "hsapiens"), i))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    save(PFM, file = c(paste0(i, ".Rda")))
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon")) != 0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
  }
}





## !!ATTENTION!! ##
# If one or more of your TFs have an abnormally large .Rda file created after the above function in /ChIPanalyser/ folder 
# then you want to check which database it gets its motif from. 
# For example, after running the above function I noticed that the TF named MAZ had an .Rda file of 24KB in size,
# much larger than any other TF which was only 1 KB.
# So this is what I did:
TFnames <- TFnames[which(TFnames == "MAZ")] # number that corresponds to MAZ
# and then I run each "if statement" separately to see when the file went down in size. I was looking for a file of 1KB.
# e.g.:
# PFM <- query(query(query(MotifDb, "hsapiens",
#                          notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10")
# save(PFM, file = c(paste0(i, ".Rda"))) 
# I checked and the file was still 24 KB.
# I tried with all datasets from the above function (the massive one that goes through all the motifs datasets)
# until I run this one from SwissRegulon that gave me the MAZ.Rda file of 1KB I was looking for. 
PFM <- query(query(query(MotifDb, "hsapiens",
                         notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon")[[1]]
save(PFM, file = c(paste0(TFnames, ".Rda")))
# If needed, manually update "PFMsource.csv" with the motif you used.
# That's it. Make sure the file was saved and move on.








## Check which TF didn't have motifs on motifDb
filename <- read.csv(paste0(cell_line,"SEbdg.txt"), sep ="\t", head = F)
TFnames <- gsub(paste0("-human_",cell_line,"_SE_trim-peaks_treat_pileup.bdg"), "", filename$V1)

PFMowned <- dir()
PFMowned <- PFMowned[grepl(".Rda", PFMowned)]
PFMowned <- gsub(".Rda", "", PFMowned)
no_motif<-TFnames[!(TFnames %in% PFMowned)]

# Here for IMR90 if I type "no_motif", I find that POLR2A doesn't have motif on motifDb. 
nomotif_tf <- "POLR2A" #change according to need
write(no_motif, file = paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/noMotif.txt"),
      append = T, sep = "\t")

# POLR2A can therefore be removed from downstream analyses.
filename <- read.csv(paste0(cell_line,"SEbdg.txt"), sep ="\t", head = F)
filename <- as.data.frame(filename[-which(filename == paste0(nomotif_tf,"-human_",cell_line,"_SE_trim-peaks_treat_pileup.bdg")),])
colnames(filename) <- "V1"
write.table(filename,
          file = paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",cell_line,"SEbdg_without",nomotif_tf,".txt"), 
          row.names = F, col.names = F, quote = F)
quit()

## Remove no motif TFs from narrowPeak datasets too
cd ..
cd ..
ls *_peaks.narrowPeak > IMR90SEnarrowPeak.txt
mv /storage/st20d/HiC_Kc167/alessandra/IMR90/preProcessing/IMR90SEnarrowPeak.txt /storage/st20d/HiC_Kc167/alessandra/IMR90/ChIPanalyser/motifs
cd ChIPanalyser/motifs/
R3.5.0

cell_line <- "IMR90" #change according to need
nomotif_tf <- "POLR2A" #change according to need
filename <- read.csv(paste0(cell_line,"SEbdg_without",nomotif_tf,".txt"), sep ="\t", head = F)
filename <- data.frame(filename[grep("treat_pileup.bdg", filename$V1),])
rownames(filename) <- NULL
TFs <- gsub(paste0("-human_",cell_line,"_SE_trim-peaks_treat_pileup.bdg"), "", filename[,1])

narrow <- read.csv(paste0(cell_line,"SEnarrowPeak.txt"),sep ='\t', head = F)
narrow <- as.data.frame(narrow[-which(narrow == paste0(nomotif_tf,"-human_",cell_line,"_SE_trim-peaks_peaks.narrowPeak")),]) # this is because POLR2A doesn't have motifs on motifDb.
colnames(narrow) <- "V1"
write.table(narrow,
            file = paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",cell_line,"SEnarrowPeak_without",nomotif_tf,".txt"), 
            row.names = F, col.names = F, quote = F)
narrow <- read.csv(paste0(cell_line,"SEnarrowPeak_without",nomotif_tf,".txt"), sep ='\t', head = F)

save(filename, narrow, TFs, file = "objectsForTable.RData")
quit()