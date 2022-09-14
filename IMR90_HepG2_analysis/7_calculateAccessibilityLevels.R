## Get accessibility percentage per each QDA for all methods of DNA accessibility ##
# In my thesis this can help me discussing why some methods (ATAC or DNase) are more valid than others (MNase)
# and that QDA=0.99 doesn't always correspond to top 1% accessibility etc.

library(rtracklayer)
cell_line <- "IMR90"
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)

#ATAC
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC/ATAC_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")

mymat1 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat1) <- quantVec
colnames(mymat1) <- "Accessibility percentage ATAC"

for(i in quantVec){
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  mysum <- sum(width(access_data))/sum(width(access))
  mymat1[rownames(mymat1)==i] <- mysum
}
save(mymat1, file="accessibilityLevels_allQDAs.Rda")


#DNase
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase/DNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")


mymat2 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat2) <- quantVec
colnames(mymat2) <- "Accessibility percentage DNase"

for(i in quantVec){
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  mysum <- sum(width(access_data))/sum(width(access))
  mymat2[rownames(mymat2)==i] <- mysum
}
save(mymat2, file="accessibilityLevels_allQDAs.Rda")



#MNase
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase/MNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")

mymat3 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat3) <- quantVec
colnames(mymat3) <- "Accessibility percentage MNase"

for(i in quantVec){
  access$score <- max(access$score) - access$score + min(access$score)
  access_data<-reduce(access[which(access$score >= quantile(access$score,(i)))])
  mysum <- sum(width(access_data))/sum(width(access))
  mymat3[rownames(mymat3)==i] <- mysum
}
save(mymat3, file="accessibilityLevels_allQDAs.Rda")



#NOMe
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/NOMe"))
get(load(file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/NOMe/NOMe_profile_",cell_line,".Rda")))

mymat4 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat4) <- quantVec
colnames(mymat4) <- "Accessibility percentage NOMe"

for(i in quantVec){
  access_data <- profile + 10 # Expand the DNA accessibility ranges - single point accessibility can be odd to work with in ChIPanalyser
  access_data<-reduce(access_data[which(access_data$score >= quantile(access_data$score,i))])
  mysum <- sum(width(access_data))/sum(width(profile + 10))
  mymat4[rownames(mymat4)==i] <- mysum
}
save(mymat4, file="accessibilityLevels_allQDAs.Rda")



#Together 
bigM <- cbind(mymat1,mymat2,mymat3,mymat4)
save(bigM, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/accessibilityLevels_allQDAs_allMethods.Rda"))





cell_line <- "HEPG2"
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)

#ATAC
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC/ATAC_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")

mymat1 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat1) <- quantVec
colnames(mymat1) <- "Accessibility percentage ATAC"

for(i in quantVec){
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  mysum <- sum(width(access_data))/sum(width(access))
  mymat1[rownames(mymat1)==i] <- mysum
}
save(mymat1, file="~/AUCs/accessibilityLevels_allQDAs.Rda")


#DNase
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase/DNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")


mymat2 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat2) <- quantVec
colnames(mymat2) <- "Accessibility percentage DNase"

for(i in quantVec){
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  mysum <- sum(width(access_data))/sum(width(access))
  mymat2[rownames(mymat2)==i] <- mysum
}
save(mymat2, file="~/AUCs/accessibilityLevels_allQDAs.Rda")



#MNase
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase"))
access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase/MNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")

mymat3 <- matrix(nrow=12,ncol=1,data=NA)
rownames(mymat3) <- quantVec
colnames(mymat3) <- "Accessibility percentage MNase"

for(i in quantVec){
  #access$score <- max(access$score) - access$score + min(access$score)
  #access_data<-reduce(access[which(access$score >= quantile(access$score,(i)))])
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase/access_",cell_line,"_",i,".Rda"))
  mysum <- sum(width(access_data))/sum(width(access))
  mymat3[rownames(mymat3)==i] <- mysum
}
save(mymat3, file="~/AUCs/accessibilityLevels_allQDAs.Rda")



#Together 
bigM <- cbind(mymat1,mymat2,mymat3)
save(bigM, file=paste0("~/AUCs/accessibilityLevels_allQDAs_allMethods_HepG2.Rda"))