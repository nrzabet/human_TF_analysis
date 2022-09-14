rm(list=ls())
library(ChIPanalyser)
setwd("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/optimalParam")
#load list of optimal QDA for each TF
optQDA <- get(load("K562OptimalQDA.Rda"))

#Get the optimal parameters for the optimal QDA for each TF
for(i in 1:nrow(optQDA)){
  setwd(c(paste0("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant", optQDA[1,])))
  ##this is just due to a typo in some of my flie names that I never got around to fixing
  ##some files are missing a "_" 
  if(file.exists(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562top10_OptimalOutputTraining.Rda")))){
  optimal <- get(load(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562top10_OptimalOutputTraining.Rda"))))
  optimalParam <- list(MSE=optimal$Optimal$OptimalMatrix$MSE,AUC= optimal$Optimal$OptimalMatrix$AUC)
}else{
  optimal <- get(load(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562_top10_OptimalOutputTraining.Rda"))))
  optimalParam <- list(MSE=optimal$Optimal$OptimalMatrix$MSE,AUC= optimal$Optimal$OptimalMatrix$AUC)
}
  #plot the heat maps for MSE and AUC
  setwd("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/optimalParam")
  pdf(file = c(paste0(rownames(optQDA[i,,drop=F]), "Heatmaps.pdf")))
  par(oma=c(0,0,3,0))
  layout(matrix(1:2,ncol=2, byrow=T),width=c(6,1.5,6,1.5),height=c(1,1))
  plotOptimalHeatMaps(optimalParam,layout=FALSE)
  dev.off()

  setwd(paste0("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant", optQDA[1,]))
  
  #get the accessibility QDA for plotting the occupancy profile
  Access <- get(load(c(paste0("/home/rtpop/ENCODE_ChIP_data/K562/DNase-seq/dataset1/PE/quantile/K562_DNase_PE_broad_q01_quant",optQDA[i,], ".Rda"))))
  
  #again, this is due to the typo in file names
  #just loading the chip profile
  if(file.exists(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562top10_ChIPTraining.Rda")))){
  eveLocusChip <- get(load(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562top10_ChIPTraining.Rda"))))
  }else{
  eveLocusChip <- get(load(c(paste0(rownames(optQDA[i,,drop=F]), "_ChIP_SE_K562top10_ChIPTraining.Rda"))))
  }
  
  #plot the occupancy profile
  setwd("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/optimalParam")
  pdf(file = c(paste0(rownames(optQDA[i,,drop=F]),"OccupancyProfile.pdf")), width =12, height = 4.5)
  plotOccupancyProfile(predictedProfile=optimal$ChIPProfiles,
                       ChIPScore=eveLocusChip,
                       chromatinState=Access,
                       occupancy=optimal$Occupancy,
                       goodnessOfFit=optimal$goodnessOfFit)
  dev.off()
}
