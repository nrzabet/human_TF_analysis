############## Pre-processing NOMe-seq data ##############
############## ############## ############## #############

# Click on links below to download them and then move to them to the cluster, IMR90/NOMe folder
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1001125&format=file&file=GSM1001125%5FIMR90%2Ehg18%2Egch%2Ewig%2Egz
# http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
# in NOMe folder: gunzip *.gz
# R3.6.0

library(rtracklayer)
cell_line <- "IMR90"

profile <- import("GSM1001125_IMR90.hg18.gch.wig")
chain <-import("hg18ToHg38.over.chain")
profile <-unlist(liftOver(profile,chain)) # this is to liftover from hg18
save(profile, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/NOMe/NOMe_profile_",cell_line,".Rda"))

quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access_data <- profile + 10 # Expand the DNA accessibility ranges - single point accessibility can be odd to work with in ChIPanalyser
  access_data <-reduce(access_data[which(access_data$score >= quantile(access_data$score,i))])
  save(access_data,file=access_cell_line)
}

 
########################## PRE-PROCESSING DONE! ########################## Now move on to ChIPanalyser.