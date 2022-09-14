### Generate QDA_IDs for MCF10 to use for ChIPanalyser


cell_line <- "MCF10"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/bdg"))

library(rtracklayer)

access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/bdg/MCF10A_CTRL_OH_treat_pileup.bdg"), format="bedGraph")
access_data<-reduce(access[which(access$score >= quantile(access$score,0))])

save(access_data,file=paste0("access_",cell_line,"_0.Rda"))