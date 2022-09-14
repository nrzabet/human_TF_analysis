rm(list=ls())
setwd("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant0.1")
library(ChIPanalyser)

##load optimal param filenames
filenames <- dir()
#getting only the validation optimal output files
filenames <- filenames[grep("MSE_OptimalOutputValidation.Rda", filenames)]
#removing the extra parts of the file name so only the TF name remains
filenames <- gsub("_ChIP_SE_K562_top10_MSE_OptimalOutputValidation.Rda", "", filenames)

#creating an empty data frame with 2 columns for TF name and AUC
df <- data.frame(TF = character(), optimalAUC = numeric(), stringsAsFactors = F)

#populating the empty data frame
for(i in filenames){
  #loading the optimal param file
  optimal <- get(load(c(paste0(i, "_ChIP_SE_K562_top10_MSE_OptimalOutputValidation.Rda"))))
  ##extracting the optimal AUC
  optimalAUC <- c(i,unname(as.list(optimal[[3]]@profiles[[1]][[1]][23])))
  df <- rbind(df, optimalAUC, stringsAsFactors = F)
}

save(df, file = "optimalAUCK562SEtop10Validationquant0.1.Rda")
