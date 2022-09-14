rm(list=ls())
setwd("C:/Users/Romana/Dropbox/Documents/K562/SE/optimalAUC")

##get filenames for each quant
files <- dir()
files <- files[which(grepl("All", files))]
data <- get(load(files[2]))
allQuantAUCV <- data[,1, drop = F]
colnames <- c("TF")
validation <- files[which(grepl("Validation", files))]

##get AUC for each QDA in one data frame
for(i in validation){
  dta <- get(load(i))
  quant <- gsub(".Rda", "", i)
  quant <- gsub("optimalAUCK562AllTF", '', quant)
  quant <- gsub("uant", "", quant)
  allQuantAUCV <- cbind(allQuantAUCV, dta[,2])
  colnames <- c(colnames, quant)
}
colnames(allQuantAUCV) <- colnames
rownames(allQuantAUCV) <- data[,1]
allQuantAUCV <- allQuantAUCV[,-1]
allQuantAUCV <- allQuantAUCV[,c(12,1:11)]

##figure out how many clusters with the elbow method
wss <- function(dat, k, nstart) {
  kmeans(dat, k, nstart = nstart )$tot.withinss
}
k.values <- 1:15
elbowV <- c()
for(i in k.values){
  set.seed(13)
  elbowV <- c(elbowV,wss(dat=allQuantAUCV, k=i, nstart=20))
}
#elbow plot
pdf(file="K562kmeansElbowPlotVal.pdf", height = 7, width = 6)
plot(k.values, elbowV,
       type="b", pch = 19, frame = FALSE,
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares",
       main = "Elbow plot")
dev.off()

#run k-means with k=4
set.seed(13)
clustersV <- kmeans(allQuantAUCV, 4, nstart = 20)

#colour blind friendly colours
cols_contrast <- cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 

#get the 4 clusters into separate variables
#note that the order of them might get swiched around with each re-run
#this might affect the plots
#clusters stay the same, but the order may be different
#usually the "traditional" and "partial" clusters swap 
#i.e. clust1 is sometimes the partial cluster, sometimes the traditional one and likewise for clust3
#occasionally the others swap around too
#so the plot labels might need to be tweaked to accommodate this
#annoying but haven't found a way around it
clust1V <- as.matrix(allQuantAUCV[which(grepl("1", clustersV[[1]])),])
clust2V <- as.matrix(allQuantAUCV[which(grepl("2", clustersV[[1]])),])
clust3V <- as.matrix(allQuantAUCV[which(grepl("3", clustersV[[1]])),])
clust4V <- as.matrix(allQuantAUCV[which(grepl("4", clustersV[[1]])),])
QDA <- c(seq(0, 0.9, by = 0.1), 0.95,0.99)

#plots
pdf(file = "TraditionalKmeansAllQDATVal60.pdf")
plot(x=QDA, y=clustersV[[2]][1,],type = "l", lty = 1, col = cols_contrast[6],
    xlab = "QDA", ylab = "optimal AUC",lwd = 3, xlim = c(0.1,1), ylim = c(0,1),
    main = "Traditional TF")
for(i in 1:nrow(clust1V)){
  lines(x=QDA, y=clust1V[i,], type = "l", lty = 1, col = cols_contrast[3], lwd = 0.75)
}
lines(x=QDA, y=clustersV[[2]][1,],type = "l", lty = 1, col = cols_contrast[6], lwd =3)
dev.off()

pdf(file = "PioneerKmeansAllQDATVal60.pdf")
plot(x=QDA, y=clustersV[[2]][2,],type = "l", lty = 1, col = cols_contrast[6],
    xlab = "QDA", ylab = "optimal AUC",lwd = 3, xlim = c(0.1,1), ylim = c(0,1),
    main = "Pioneer TF")
for(i in 1:nrow(clust2V)){
  lines(x=QDA, y=clust2V[i,], type = "l", lty = 1, col = cols_contrast[3], lwd = 0.75)
}
lines(x=QDA, y=clustersV[[2]][2,],type = "l", lty = 1, col = cols_contrast[6], lwd =3)
dev.off()

pdf(file = "PartialKmeansAllQDATVal60.pdf")
plot(x=QDA, y=clustersV[[2]][3,],type = "l", lty = 1, col = cols_contrast[6],
    xlab = "QDA", ylab = "optimal AUC",lwd = 3, xlim = c(0.1,1), ylim = c(0,1),
    main = "Partial pioneer TF")
for(i in 1:nrow(clust3V)){
  lines(x=QDA, y=clust3V[i,], type = "l", lty = 1, col = cols_contrast[3], lwd = 0.75)
}
lines(x=QDA, y=clustersV[[2]][3,],type = "l", lty = 1, col = cols_contrast[6], lwd =3)
dev.off()

pdf(file = "BadKmeansAllQDATVal60.pdf")
plot(x=QDA, y=clustersV[[2]][4,],type = "l", lty = 1, col = cols_contrast[6],
    xlab = "QDA", ylab = "optimal AUC",lwd = 3, xlim = c(0.1,1), ylim = c(0,1),
  main = "Poorly predicted TF")
for(i in 1:nrow(clust4V)){
  lines(x=QDA, y=clust4V[i,], type = "l", lty = 1, col = cols_contrast[3], lwd = 0.75)
}
lines(x=QDA, y=clustersV[[2]][4,],type = "l", lty = 1, col = cols_contrast[6], lwd =3)
dev.off()

#barplot of how many TFs in each cluster
#again, this suffers from the cluster switcharoo
#so make sure not to mislabel
oneV <- clustersV$cluster[grep("1",clustersV$cluster)]
twoV <- clustersV$cluster[grep("2",clustersV$cluster)]
threeV <- clustersV$cluster[grep("3",clustersV$cluster)]
fourV <- clustersV$cluster[grep("4",clustersV$cluster)]
numberV <- c(length(oneV), length(twoV), length(threeV), length(fourV))
numberV <- numberV[order(numberV, decreasing = T)]

pdf("K562SETFsinClustAllQDATVal60.pdf")
barplot(numberV, names.arg = c("Partial pioneer","Pioneer","Traditional", "Poorly predicted"), col = cols_contrast[6],
        ylim=c(0,70), ylab = "No. of TFs", main = "TF cluster sizes")
dev.off()
