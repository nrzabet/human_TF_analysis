rm(list=ls())
setwd("/home/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant0.99")
library(ChIPanalyser)

##load optimal param filenames
filenames <- dir()
#getting only the training optimal parameters files
filenames <- filenames[grep("OptimalOutputTraining.Rda", filenames)]
#removing the extra parts of the file name so only the TF name remains
filenames <- gsub("_ChIP_SE_K562top10_OptimalOutputTraining.Rda", "", filenames)

#creating an empty data frame
df <- data.frame(TF = character(), lambda = numeric(), BM = numeric(), stringsAsFactors = F)

#populating the empty data frame
for(i in filenames){
  #loading the optimal param file
  optimal <- get(load(c(paste0(i, "_ChIP_SE_K562top10_OptimalOutputTraining.Rda"))))
  ##getting the optimal parameters
  optimalParam <- c(i, unname(as.list(optimal[[1]][[1]][[12]])))
  df <- rbind(df, optimalParam, stringsAsFactors = F)
}
colnames(df) <- c("TF", "lambda", "BM")

save(df, file = "optimalParamK562PEtop10quant0.99.Rda")
