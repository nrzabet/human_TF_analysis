
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("~/ChIPanalyser_human/romana/")

################################################################################
# preprocessing -K-means - K562
################################################################################

ChIP_data_stats <- read.table("data/Table_S2.tsv", header=T, sep="\t")
ChIP_data_stats <- read.csv("data/Table_A2.csv")

ChIP_data_stats_core <- ChIP_data_stats[ChIP_data_stats$JASPAR.motif=="Y",]
ChIP_data_stats_others <- ChIP_data_stats[ChIP_data_stats$JASPAR.motif=="N" & ChIP_data_stats$Other=="Y",]
ChIP_data_stats_not <- ChIP_data_stats[ChIP_data_stats$JASPAR.motif=="N" & ChIP_data_stats$Other=="N",]

write.table(ChIP_data_stats_core, file="data/Table_S2.csv", sep=",", row.names = FALSE)

pdf("Figure_S1_v002.pdf", width=12, height=6,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2,ncol=1,byrow=T), width = c(nrow(ChIP_data_stats_core)+5),height = c(1,1))

par(mar=c(4, 5.0, 3.5, 0)+0.1);

ord <- order(ChIP_data_stats_core$Overall.alignment, decreasing = TRUE)
m <- barplot(ChIP_data_stats_core$Overall.alignment[ord], main="Overall Alignment (%) - JASPAR core", xlab="", ylab="", col=cbbPalette[3], xaxt="n", ylim=c(0,100))
axis(1, at=m, labels=ChIP_data_stats_core$TF.name[ord], col.axis="black", las=2, cex.axis=0.5)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)

ord <- order(ChIP_data_stats_core$No..of.peaks, decreasing = TRUE)
m <- barplot(log10(ChIP_data_stats_core$No..of.peaks[ord]), main="Number of peaks - JASPAR core", xlab="", ylab="", col=cbbPalette[3], xaxt="n", yaxt="n", ylim=c(0,6))
axis(1, at=m, labels=ChIP_data_stats_core$TF.name[ord], col.axis="black", las=2, cex.axis=0.5)
axis(2, at=0:6, labels=10^(0:6), col.axis="black", las=2, cex.axis=1.0)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
dev.off()




################################################################################
# QDA analysis -K-means - K562
################################################################################
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
QDA_ID <- 1:12

# Load Romana's data
# load("AUCs/K562AllQuantAUCAllTFVal_romana.RData")
# allQuantAUCV <- as.matrix(allQuantAUCV)
AUCmatrix_K562 <- read.table("data/Table_S3.tsv", header=T, sep="\t",row.names = 1)
colnames(AUCmatrix_K562) <- QDA_ID

AUCmatrix_K562_core <- AUCmatrix_K562[rownames(AUCmatrix_K562) %in% ChIP_data_stats_core$TF.name,]

write.table(AUCmatrix_K562_core, file="data/Table_S3.csv", sep=",", row.names = TRUE)


rownames(AUCmatrix_K562_core) <- paste0(rownames(AUCmatrix_K562_core),"_K562")


wss <- function(dat, k, nstart){
  kmeans(dat, k, nstart = nstart)$tot.withinss
}

k.values <- 1:15 # why 15?
elbowV <- c()

for(i in k.values){
  set.seed(13)
  elbowV <- c(elbowV,wss(dat=AUCmatrix_K562_core[,1:11], k=i, nstart=20)) # why start at 20?
}

# Elbow plot
pdf(file=paste0("KMeans_ElbowPlot_K562_core.pdf"), height=5, width=6)
plot(k.values, elbowV,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main = "Elbow plot",
     ylim=c(0,20),
     xlim=c(0,15))
dev.off()

colfunc<-colorRampPalette(c("white",cbbPalette[6]))
Colors <- colfunc(20)

#
# k <- "kValue3"
# set.seed(13)
# clustersV <- kmeans(AUCmatrix_K562_core[,1:11], 3, nstart=20)
#
# clust1V <- as.matrix(AUCmatrix_K562_core[which(grepl("1", clustersV[[1]])),])
# clust2V <- as.matrix(AUCmatrix_K562_core[which(grepl("2", clustersV[[1]])),])
# clust3V <- as.matrix(AUCmatrix_K562_core[which(grepl("3", clustersV[[1]])),])
#
# clust1V_K562_core <- clust1V
# clust2V_K562_core <- clust2V
# clust3V_K562_core <- clust3V
#
# rownames(clust1V_K562_core) <- gsub("_K562", "", rownames(clust1V_K562_core))
# rownames(clust2V_K562_core) <- gsub("_K562", "", rownames(clust2V_K562_core))
# rownames(clust3V_K562_core) <- gsub("_K562", "", rownames(clust3V_K562_core))
#
# colfunc<-colorRampPalette(c("white",cbbPalette[6]))
# Colors <- colfunc(20)
#
#
# # all
# pdf(file = paste0("Kmeans_K562_core_clusters_heatmap.pdf"), width=6, height=20,pointsize = 12)
# par(cex=1.1);
# par(las=1)
# options(scipen=999)
#
# layout(matrix(1:3, ncol=1),
#        height=c(nrow(clust1V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)),
#                 nrow(clust2V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)),
#                 nrow(clust3V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)) ) )
# par(mar=c(4, 9.0, 3.5, 0)+0.1);
#
# # cluster 1
# clust <- clust1V_K562_core[order(rownames(clust1V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("ADFs (",nrow(clust1V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)
#
#
#
# # cluster 2
# clust <- clust2V_K562_core[order(rownames(clust2V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("partial AIFs (",nrow(clust2V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
#
#
# # cluster 3
# clust <- clust3V_K562_core[order(rownames(clust3V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("AIFs (",nrow(clust3V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)
# dev.off()
#
#
#
#
# pdf(file = paste0("Kmeans_K562_3clusters_heatmap_t_legend.pdf"), width=0.75, height=5,pointsize = 12)
# legend_image<-as.raster(matrix(rev(Colors),ncol=1))
# # raster scacle
# par(mar=c(3,0.5,3.2,0.5))
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
# text(x=1.6, y =seq(0,1,l=6) , labels = round(seq(0.5,1.0,l=6),2),cex=0.95)
# rasterImage(legend_image, 0, 0, 1,1)
# dev.off()
#
#

pdf("Figure_S2_maxAUC_K562_core.pdf", width=5, height=4,pointsize = 12)
m <- hist(apply(AUCmatrix_K562_core,1,max ), col=cbbPalette[3], breaks=seq(0,1, by=0.1), main="Maximum AUC", ylab="count", xlab="AUC", ylim=c(0,70))
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)
dev.off()

sum(apply(AUCmatrix_K562_core,1,max ) > 0.8)/nrow(AUCmatrix_K562_core)





k <- "kValue2"
set.seed(13)
clustersV <- kmeans(AUCmatrix_K562_core[,1:11], 2, nstart=20)

clust1V <- as.matrix(AUCmatrix_K562_core[which(grepl("1", clustersV[[1]])),])
clust2V <- as.matrix(AUCmatrix_K562_core[which(grepl("2", clustersV[[1]])),])

clust1V_K562_core <- clust1V
clust2V_K562_core <- clust2V

rownames(clust1V_K562_core) <- gsub("_K562", "", rownames(clust1V_K562_core))
rownames(clust2V_K562_core) <- gsub("_K562", "", rownames(clust2V_K562_core))


# all
pdf(file = paste0("Kmeans_K562_core_2clusters_heatmap.pdf"), width=9, height=23,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2, ncol=1),
       height=c(nrow(clust1V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core)),
                nrow(clust2V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core))) )
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 1
clust <- clust1V_K562_core[order(rownames(clust1V_K562_core)),]
image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("ADFs (",nrow(clust1V_K562_core),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)



# cluster 2
clust <- clust2V_K562_core[order(rownames(clust2V_K562_core)),]
image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs (",nrow(clust2V_K562_core),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
dev.off()




################################################################################
# QDA analysis - threshold - K562
################################################################################


#
TF_high <- apply(AUCmatrix_K562_core[,9:11],1,mean, na.rm=TRUE)
TF_low <- apply(AUCmatrix_K562_core[,1:3],1,mean, na.rm=TRUE)
TF_diff <- TF_high - TF_low
manual_class <- rep("other", nrow(AUCmatrix_K562_core))
manual_class[which((TF_low >= 0.8 | TF_high >= 0.8) & abs(TF_diff) < 0.1)] <- "AIF"
manual_class[which((TF_high >= 0.8) & TF_diff >= 0.3)] <- "ADF"
manual_class[which((TF_low >= 0.8) & TF_diff <= -0.3)] <- "IDF"
manual_class[which(!(TF_low >= 0.65 | TF_high >= 0.65))] <- "poorly predicted"
manual_class[which((TF_high >= 0.65) & (TF_diff >= 0.1 & TF_diff < 0.3))] <- "partial AIF/ADF"
manual_class[which((TF_low >= 0.65) & (TF_diff < -0.1 & TF_diff >= -0.3))] <- "partial AIF/IDF"
table(manual_class)
split(AUCmatrix_K562_core, manual_class)
#
# TF_high <- apply(AUCmatrix_K562_core[,10:11],1,mean, na.rm=TRUE)
# TF_low <- apply(AUCmatrix_K562_core[,1:2],1,mean, na.rm=TRUE)
# TF_diff <- TF_high - TF_low
# manual_class <- rep("none", nrow(AUCmatrix_K562_core))
# manual_class[which(TF_low < 0.8 & TF_high >= 0.8 & TF_diff < 0.2)] <- "partial AIF/ADF"
# manual_class[which(TF_low < 0.8 & TF_high >= 0.8 & TF_diff >= 0.2)] <- "ADF"
# manual_class[which(TF_low >= 0.8 & TF_high < 0.8 & abs(TF_diff < 0.2))] <- "partial AIF/IDF"
# manual_class[which(TF_low >= 0.8 & TF_high >= 0.8 & TF_diff < 0.2)] <- "AIF"
# manual_class[which(TF_low >= 0.8 & TF_high < 0.8 & abs(TF_diff >= 0.2))] <- "IDF"
# table(manual_class)


rownames(AUCmatrix_K562_core) <- gsub("_K562", "", rownames(AUCmatrix_K562_core))
AUCmatrix_K562_core_manual_class <- split(AUCmatrix_K562_core[order(TF_high+TF_diff),], manual_class[order(TF_high+TF_diff)])

# all
pos <- 1:12
pdf(file = paste0("Kmeans_K562_core_threshold_heatmap.pdf"), width=9, height=30,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:6, ncol=1),
       height=c(nrow(AUCmatrix_K562_core_manual_class[["ADF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["AIF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["other"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["poorly predicted"]])/nrow(AUCmatrix_K562_core)))
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 1
#clust <- AUCmatrix_K562_core_manual_class[["ADF"]][order(rownames(AUCmatrix_K562_core_manual_class[["ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["ADF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)



# cluster 2
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][order(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)


# cluster 3
#clust <- AUCmatrix_K562_core_manual_class[["AIF"]][order(rownames(AUCmatrix_K562_core_manual_class[["AIF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["AIF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["AIF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["AIF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("AIFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)

# cluster 4
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][order(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/IDFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[4], line = 0.7, adj = -0.09, cex=1.6)

# cluster 5
#clust <- AUCmatrix_K562_core_manual_class[["other"]][order(rownames(AUCmatrix_K562_core_manual_class[["other"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["other"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["other"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["other"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("others (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[5], line = 0.7, adj = -0.09, cex=1.6)

# cluster 6
#clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][order(rownames(AUCmatrix_K562_core_manual_class[["poorly predicted"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["poorly predicted"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("poorly predicted (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[6], line = 0.7, adj = -0.09, cex=1.6)

dev.off()





# all
pos <- 1:12
pdf(file = paste0("Kmeans_K562_core_threshold_heatmap_large.pdf"), width=9, height=25,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2, ncol=1),
       height=c(nrow(AUCmatrix_K562_core_manual_class[["AIF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]])/nrow(AUCmatrix_K562_core)))
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 3
#clust <- AUCmatrix_K562_core_manual_class[["AIF"]][order(rownames(AUCmatrix_K562_core_manual_class[["AIF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["AIF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["AIF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["AIF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("AIFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)

# cluster 2
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][order(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)

dev.off()



pdf(file = paste0("Kmeans_K562_core_threshold_heatmap_small.pdf"), width=9, height=12.5,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:4, ncol=1),
       height=c(nrow(AUCmatrix_K562_core_manual_class[["ADF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["other"]])/nrow(AUCmatrix_K562_core),
                nrow(AUCmatrix_K562_core_manual_class[["poorly predicted"]])/nrow(AUCmatrix_K562_core)))
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 1
#clust <- AUCmatrix_K562_core_manual_class[["ADF"]][order(rownames(AUCmatrix_K562_core_manual_class[["ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["ADF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)



# cluster 4
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][order(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/IDFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[4], line = 0.7, adj = -0.09, cex=1.6)

# cluster 5
#clust <- AUCmatrix_K562_core_manual_class[["other"]][order(rownames(AUCmatrix_K562_core_manual_class[["other"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["other"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["other"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["other"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("others (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[5], line = 0.7, adj = -0.09, cex=1.6)

# cluster 6
#clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][order(rownames(AUCmatrix_K562_core_manual_class[["poorly predicted"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core_manual_class[["poorly predicted"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][data.order,pos]
clust <- AUCmatrix_K562_core_manual_class[["poorly predicted"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("poorly predicted (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[6], line = 0.7, adj = -0.09, cex=1.6)

dev.off()























#
# # all
# pos <- 2:10
# pdf(file = paste0("Kmeans_K562_core_threshold_heatmap.pdf"), width=9, height=23,pointsize = 12)
# par(cex=1.1);
# par(las=1)
# options(scipen=999)
# # cluster 3
# #clust <- AUCmatrix_K562_core_manual_class[["AIF"]][order(rownames(AUCmatrix_K562_core_manual_class[["AIF"]])),pos]
# data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_K562_core[,pos])))
# data.order <- order.dendrogram(data.dendro)
# clust <- AUCmatrix_K562_core[data.order,pos]
#
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("AIFs (",nrow(clust),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)
#
#
#
# dev.off()
#
#
#
#
#
# wss <- function(dat, k, nstart){
#   kmeans(dat, k, nstart = nstart)$tot.withinss
# }
#
# k.values <- 1:15 # why 15?
# elbowV <- c()
#
# for(i in k.values){
#   set.seed(13)
#   elbowV <- c(elbowV,wss(dat=AUCmatrix_K562_core[,2:11], k=i, nstart=20)) # why start at 20?
# }
#
# # Elbow plot
# pdf(file=paste0("KMeans_ElbowPlot_K562_core_withoutextremes.pdf"), height=5, width=6)
# plot(k.values, elbowV,
#      type="b", pch = 19, frame = FALSE,
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares",
#      main = "Elbow plot",
#      ylim=c(0,20),
#      xlim=c(0,15))
# dev.off()
#
#
#
# k <- "kValue2"
# set.seed(13)
# clustersV <- kmeans(AUCmatrix_K562_core[,pos], 3, nstart=20)
#
# clust1V <- as.matrix(AUCmatrix_K562_core[which(grepl("1", clustersV[[1]])),pos])
# clust2V <- as.matrix(AUCmatrix_K562_core[which(grepl("2", clustersV[[1]])),pos])
# clust3V <- as.matrix(AUCmatrix_K562_core[which(grepl("3", clustersV[[1]])),pos])
#
#
# clust1V_K562_core <- clust1V
# clust2V_K562_core <- clust2V
# clust3V_K562_core <- clust3V
#
# rownames(clust1V_K562_core) <- gsub("_K562", "", rownames(clust1V_K562_core))
# rownames(clust2V_K562_core) <- gsub("_K562", "", rownames(clust2V_K562_core))
# rownames(clust3V_K562_core) <- gsub("_K562", "", rownames(clust3V_K562_core))
#
# colfunc<-colorRampPalette(c("white",cbbPalette[6]))
# Colors <- colfunc(20)
#
#
# # all
# pdf(file = paste0("Kmeans_K562_core_clusters_heatmap_withoutextremes.pdf"), width=6, height=20,pointsize = 12)
# par(cex=1.1);
# par(las=1)
# options(scipen=999)
#
# layout(matrix(1:3, ncol=1),
#        height=c(nrow(clust1V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)),
#                 nrow(clust2V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)),
#                 nrow(clust3V_K562_core)/nrow(rbind(clust1V_K562_core, clust2V_K562_core, clust3V_K562_core)) ) )
# par(mar=c(4, 9.0, 3.5, 0)+0.1);
#
# # cluster 1
# clust <- clust1V_K562_core[order(rownames(clust1V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("ADFs (",nrow(clust1V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)
#
#
#
# # cluster 2
# clust <- clust2V_K562_core[order(rownames(clust2V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("partial AIFs (",nrow(clust2V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
#
#
# # cluster 3
# clust <- clust3V_K562_core[order(rownames(clust3V_K562_core)),]
# image(1:ncol(clust),1:nrow(clust),t(clust),
#       axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
# title(main=paste0("AIFs (",nrow(clust3V_K562_core),")"),cex.main=1.8)
# #title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
# title(xlab="QDAs",line=4.5, cex.lab=1.2)
# axis(1, at=QDA_ID, labels=quantVec,las = 1,cex.axis=1.2)
# axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
# mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)
# dev.off()



################################################################################
# preprocessing - mm10
################################################################################

AUCmatrix_mm10 <- read.csv("data/Table_S5.csv", header=T, row.names = 1)
colnames(AUCmatrix_mm10) <- QDA_ID
library(MotifDb)

motifs <- query(MotifDb, andStrings=c("CTCF", "musculus"),
                orStrings=c("jaspar"),
                notStrings="ctcfl")

CTCF.indices = grep ('^CTCF', values (MotifDb)$geneSymbol, ignore.case=TRUE)

jaspar.CTCF.matrices <- subset (MotifDb, geneSymbol=='CTCF' &
                                  dataSource == 'JASPAR_CORE')

rowTFs <- unlist(strsplit(rownames(AUCmatrix_mm10), "_"))[c(T,F)]
rowCell <- unlist(strsplit(rownames(AUCmatrix_mm10), "_"))[c(F,T)]
mm10_preprocess_stats <- data.frame("TF"=rowTFs, "cell"=rowCell, "alignement"=0, "peaks"=0)

jaspar_core_hits <- function(x){
  x <- query(MotifDb, andStrings=c(x, "jaspar"), orStrings=c("Mmusculus", "Hsapiens"))
  return(length(x)>0)
}


jaspar_motif <- sapply(rowTFs, jaspar_core_hits)
corrected <- c("MAZ", "TBP", "ZKSCAN1")
jaspar_motif[which(names(jaspar_motif)%in%corrected)] <- TRUE
to_remove <- c(which(((mm10_preprocess_stats$TF=="GATA1" & mm10_preprocess_stats$cell=="G1E") |
                        (mm10_preprocess_stats$TF=="MYOG" & mm10_preprocess_stats$cell=="C2C12"))))
jaspar_motif[to_remove] <- FALSE
AUCmatrix_mm10_core <- AUCmatrix_mm10[jaspar_motif,]
write.table(AUCmatrix_mm10_core, file="data/Table_S7_core.csv", sep=",", row.names = TRUE)

mm10_preprocess_stats_core <- mm10_preprocess_stats[jaspar_motif,]

mm10_align_all <- read.csv("data/mm10_align_all.csv", header=T, row.names = 1)
mm10_peaks_all <- read.csv("data/mm10_peaks_all.csv", header=T, row.names = 1)

for(i in 1:nrow(mm10_preprocess_stats_core)){
  row_id <- which(rownames(mm10_align_all)==mm10_preprocess_stats_core[i,1])
  col_id <- which(colnames(mm10_align_all)==gsub("-",".",mm10_preprocess_stats_core[i,2]))
  mm10_preprocess_stats_core[i, 3] <- mm10_align_all[row_id, col_id]

  row_id <- which(rownames(mm10_peaks_all)==mm10_preprocess_stats_core[i,1])
  col_id <- which(colnames(mm10_peaks_all)==gsub("-",".",mm10_preprocess_stats_core[i,2]))
  mm10_preprocess_stats_core[i, 4] <- mm10_peaks_all[row_id, col_id]
}

write.table(mm10_preprocess_stats_core, file="data/Table_S5_core.csv", sep=",", row.names = TRUE)


cellsFactors <-as.factor(mm10_preprocess_stats_core$cell)

pdf("Figure_S4.pdf", width=12, height=6,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2,ncol=1,byrow=T), width = c(nrow(mm10_preprocess_stats_core)+5),height = c(1,1))

par(mar=c(4, 5.0, 3.5, 0)+0.1);

cols <- cbPalette[1:8]

ord <- order(mm10_preprocess_stats_core$alignement, decreasing = TRUE)
m <- barplot(100*mm10_preprocess_stats_core$alignement[ord],
             main="Overall Alignment (%) - JASPAR core", xlab="", ylab="",
             col=cols[cellsFactors[ord]], xaxt="n", ylim=c(0,100))
axis(1, at=m, labels=paste0(mm10_preprocess_stats_core$TF[ord]," ", mm10_preprocess_stats_core$cell[ord]), col.axis="black", las=2, cex.axis=0.5)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)

ord <- order(mm10_preprocess_stats_core$peaks, decreasing = TRUE)
m <- barplot(log10(mm10_preprocess_stats_core$peaks[ord]),
             main="Number of peaks - JASPAR core", xlab="", ylab="",
             col=cols[cellsFactors[ord]], xaxt="n", yaxt="n", ylim=c(0,6))
axis(1, at=m, labels=paste0(mm10_preprocess_stats_core$TF[ord]," ", mm10_preprocess_stats_core$cell[ord]), col.axis="black", las=2, cex.axis=0.5)
axis(2, at=0:6, labels=10^(0:6), col.axis="black", las=2, cex.axis=1.0)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
legend("topright", fill=cols[1:length(levels(cellsFactors))],
       legend=as.character(levels(cellsFactors)), bty="n", cex=0.65, horiz = F)
dev.off()

################################################################################
# mm10 grouping - threshold
################################################################################
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
QDA_ID <- 1:12

#
TF_high <- apply(AUCmatrix_mm10_core[,9:11],1,mean, na.rm=TRUE)
TF_low <- apply(AUCmatrix_mm10_core[,1:3],1,mean, na.rm=TRUE)
TF_diff <- TF_high - TF_low
manual_class <- rep("other", nrow(AUCmatrix_mm10_core))
manual_class[which((TF_low >= 0.8 | TF_high >= 0.8) & abs(TF_diff) < 0.1)] <- "AIF"
manual_class[which((TF_high >= 0.8) & TF_diff >= 0.3)] <- "ADF"
manual_class[which((TF_low >= 0.8) & TF_diff <= -0.3)] <- "IDF"
manual_class[which(!(TF_low >= 0.65 | TF_high >= 0.65))] <- "poorly predicted"
manual_class[which((TF_high >= 0.65) & (TF_diff >= 0.1 & TF_diff < 0.3))] <- "partial AIF/ADF"
manual_class[which((TF_low >= 0.65) & (TF_diff < -0.1 & TF_diff >= -0.3))] <- "partial AIF/IDF"
table(manual_class)

AUCmatrix_mm10_core_manual_class <- AUCmatrix_mm10_core
rownames(AUCmatrix_mm10_core_manual_class) <- gsub("_", " ", rownames(AUCmatrix_mm10_core_manual_class))
AUCmatrix_mm10_core_manual_class <- split(AUCmatrix_mm10_core_manual_class[order(TF_high+TF_diff),], manual_class[order(TF_high+TF_diff)])

colfunc<-colorRampPalette(c("white",cbbPalette[6]))
Colors <- colfunc(20)

# all
pos <- 1:12
pdf(file = paste0("Kmeans_mm10_core_threshold_heatmap_large.pdf"), width=9, height=15,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2, ncol=1),
       height=c(nrow(AUCmatrix_mm10_core_manual_class[["AIF"]])/nrow(AUCmatrix_mm10_core),
                nrow(AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]])/nrow(AUCmatrix_mm10_core)))
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 3
#clust <- AUCmatrix_mm10_core_manual_class[["AIF"]][order(rownames(AUCmatrix_mm10_core_manual_class[["AIF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["AIF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["AIF"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["AIF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("AIFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)

# cluster 2
#clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]][order(rownames(AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)

dev.off()



pdf(file = paste0("Kmeans_mm10_core_threshold_heatmap_small.pdf"), width=9, height=10,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:4, ncol=1),
       height=c(nrow(AUCmatrix_mm10_core_manual_class[["ADF"]])*0.7/nrow(AUCmatrix_mm10_core),
                nrow(AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]])/nrow(AUCmatrix_mm10_core),
                nrow(AUCmatrix_mm10_core_manual_class[["other"]])/nrow(AUCmatrix_mm10_core),
                nrow(AUCmatrix_mm10_core_manual_class[["poorly predicted"]])/nrow(AUCmatrix_mm10_core)))
par(mar=c(4, 9.0, 3.5, 0)+0.1);

# cluster 1
#clust <- AUCmatrix_mm10_core_manual_class[["ADF"]][order(rownames(AUCmatrix_mm10_core_manual_class[["ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["ADF"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)



# cluster 4
#clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]][order(rownames(AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/IDFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[4], line = 0.7, adj = -0.09, cex=1.6)

# cluster 5
#clust <- AUCmatrix_mm10_core_manual_class[["other"]][order(rownames(AUCmatrix_mm10_core_manual_class[["other"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["other"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["other"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["other"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("others (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[5], line = 0.7, adj = -0.09, cex=1.6)

# cluster 6
#clust <- AUCmatrix_mm10_core_manual_class[["poorly predicted"]][order(rownames(AUCmatrix_mm10_core_manual_class[["poorly predicted"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_mm10_core_manual_class[["poorly predicted"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_mm10_core_manual_class[["poorly predicted"]][data.order,pos]
clust <- AUCmatrix_mm10_core_manual_class[["poorly predicted"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("poorly predicted (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[6], line = 0.7, adj = -0.09, cex=1.6)

dev.off()


pdf("Figure_4_maxAUC_mm10_core.pdf", width=5, height=4,pointsize = 12)
m <- hist(apply(AUCmatrix_mm10_core,1,max ), col=cbbPalette[3],
          breaks=seq(0,1, by=0.1), main="Maximum AUC", ylab="count", xlab="AUC", ylim=c(0,50))
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)
dev.off()

sum(apply(AUCmatrix_mm10_core,1,max ) > 0.8)/nrow(AUCmatrix_mm10_core)


################################################################################
# mm10 counts - threshold
################################################################################
mm10_counts <- table(mm10_preprocess_stats_core$TF)

K562_names <- ChIP_data_stats_core$TF.name
K562_names <- gsub("3xFLAG-","",K562_names)
K562_names <- gsub("eGFP-","",K562_names)
K562_counts <- table(K562_names)



ids <- match(names(mm10_counts),names(K562_counts))
mm10_K562_counts <- mm10_counts
mm10_K562_counts[!is.na(ids)] <- mm10_K562_counts[!is.na(ids)] + K562_counts[ids[!is.na(ids)]]

mm10_K562_counts <- mm10_K562_counts[order(mm10_K562_counts, decreasing = T)]

mm10_K562_counts_manual <- matrix(0, ncol=length(mm10_K562_counts), nrow=6)
rownames(mm10_K562_counts_manual) <- c("AIFs", "partial AIFs/ADFs", "ADFs", "partial AIFs/IDFs", "others", "poorly predicted")
colnames(mm10_K562_counts_manual) <- names(mm10_K562_counts)

for(i in 1:ncol(mm10_K562_counts_manual)){
  TF <- colnames(mm10_K562_counts_manual)[i]
  mm10_K562_counts_manual[1,i] <- mm10_K562_counts_manual[1,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["AIF"]])))
  mm10_K562_counts_manual[2,i] <- mm10_K562_counts_manual[2,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["partial AIF/ADF"]])))
  mm10_K562_counts_manual[3,i] <- mm10_K562_counts_manual[3,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["ADF"]])))
  mm10_K562_counts_manual[4,i] <- mm10_K562_counts_manual[4,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["partial AIF/IDF"]])))
  mm10_K562_counts_manual[5,i] <- mm10_K562_counts_manual[5,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["other"]])))
  mm10_K562_counts_manual[6,i] <- mm10_K562_counts_manual[6,i] + length(grep(paste0(TF," "), rownames(AUCmatrix_mm10_core_manual_class[["poorly predicted"]])))

  TF_variations <- c(TF, paste0("3xFLAG-", TF), paste0("eGFP-",TF))

  mm10_K562_counts_manual[1,i] <- mm10_K562_counts_manual[1,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["AIF"]]) %in% TF_variations)
  mm10_K562_counts_manual[2,i] <- mm10_K562_counts_manual[2,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/ADF"]]) %in% TF_variations)
  mm10_K562_counts_manual[3,i] <- mm10_K562_counts_manual[3,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["ADF"]]) %in% TF_variations)
  mm10_K562_counts_manual[4,i] <- mm10_K562_counts_manual[4,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["partial AIF/IDF"]]) %in% TF_variations)
  mm10_K562_counts_manual[5,i] <- mm10_K562_counts_manual[5,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["other"]]) %in% TF_variations)
  mm10_K562_counts_manual[6,i] <- mm10_K562_counts_manual[6,i] + sum(rownames(AUCmatrix_K562_core_manual_class[["poorly predicted"]]) %in% TF_variations)
}

pdf("Figure_4_classifications.pdf", width=8, height=4,pointsize = 12)
cols <- c(cbPalette[c(6,4,3,8,1,2)])
#cols <- c(cbPalette[c(7,5,4,2,3,1)])
cols <- c(cbPalette[c(6,5,4,8,2,1)])
cols <- c(cbPalette[c(6,4,5,8,2,1)])
par(las=2)
barplot(mm10_K562_counts_manual, col=cols, main="classifications of TFs in mouse cell lines and human K562 cells", ylab="datasets")
legend("topright", fill=cols, legend = rownames(mm10_K562_counts_manual), bty="n")
dev.off()



################################################################################
# preprocessing - IMR90
################################################################################

IMR90_HepG2_preprocess_stats <- read.csv("data/Table_S8_all.csv", header=T)

rowTFs <- IMR90_HepG2_preprocess_stats[,1]
rowCell <- IMR90_HepG2_preprocess_stats[,2]
#
# mm10_preprocess_stats <- data.frame("TF"=rowTFs, "cell"=rowCell, "alignement"=0, "peaks"=0)
# AUCmatrix_mm10_core <- AUCmatrix_mm10[jaspar_motif,]
# write.table(AUCmatrix_mm10_core, file="data/Table_S7_core.csv", sep=",", row.names = TRUE)

IMR90_HepG2_preprocess_stats_core <- IMR90_HepG2_preprocess_stats[IMR90_preprocess_stats$JASPAR=="Y",]

#write.table(IMR90_HepG2_preprocess_stats_core, file="data/Table_S8_core.csv", sep=",", row.names = FALSE)


IMR90_HepG2_preprocess_stats_core <- read.csv("data/Table_S8_core.csv", header=T)

cellsFactors <- as.factor(IMR90_HepG2_preprocess_stats_core$Cell.line)

cols <- cbPalette[c(2,3)]

pdf("Figure_S5_TFs.pdf", width=4.7, height=6,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2,ncol=1,byrow=T), width = c(nrow(mm10_preprocess_stats_core)+5),height = c(1,1))

par(mar=c(4, 5.0, 3.5, 0)+0.1);


ord <- order(IMR90_HepG2_preprocess_stats_core$Percentage, decreasing = TRUE)
m <- barplot(IMR90_HepG2_preprocess_stats_core$Percentage[ord],
             main="Overall Alignment (%) - JASPAR core", xlab="", ylab="",
             col=cols[cellsFactors[ord]], xaxt="n", ylim=c(0,100))
axis(1, at=m, labels=paste0(IMR90_HepG2_preprocess_stats_core$TFs[ord]," ", IMR90_HepG2_preprocess_stats_core$Cell.line[ord]), col.axis="black", las=2, cex.axis=0.5)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)

ord <- order(IMR90_HepG2_preprocess_stats_core$ChIP.seq.peaks.number, decreasing = TRUE)
m <- barplot(log10(IMR90_HepG2_preprocess_stats_core$ChIP.seq.peaks.number[ord]),
             main="Number of peaks - JASPAR core", xlab="", ylab="",
             col=cols[cellsFactors[ord]], xaxt="n", yaxt="n", ylim=c(0,6))
axis(1, at=m, labels=paste0(IMR90_HepG2_preprocess_stats_core$TFs[ord]," ", IMR90_HepG2_preprocess_stats_core$Cell.line[ord]), col.axis="black", las=2, cex.axis=0.5)
axis(2, at=0:6, labels=10^(0:6), col.axis="black", las=2, cex.axis=1.0)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)
legend("topright", fill=cols[1:length(levels(cellsFactors))],
       legend=as.character(levels(cellsFactors)), bty="n", cex=0.65, horiz = F)
dev.off()


load("objects_Ale/AUCs_allRegions_IMR90_DNase.RData")
rownames(allAUCs) <- paste0(rownames(allAUCs), "_DNase_IMR90")
AUCmatrix_IMR90 <- allAUCs
load("objects_Ale/AUCs_allRegions_IMR90_ATAC.RData")
rownames(allAUCs) <- paste0(rownames(allAUCs), "_ATAC_IMR90")
AUCmatrix_IMR90 <- rbind(AUCmatrix_IMR90, allAUCs)
load("objects_Ale/AUCs_allRegions_IMR90_MNase.RData")
rownames(allAUCs) <- paste0(rownames(allAUCs), "_MNase_IMR90")
AUCmatrix_IMR90 <- rbind(AUCmatrix_IMR90, allAUCs)
load("objects_Ale/AUCs_allRegions_IMR90_NOMe.RData")
rownames(allAUCs) <- paste0(rownames(allAUCs), "_NOMe_IMR90")
AUCmatrix_IMR90 <- rbind(AUCmatrix_IMR90, allAUCs)

AUCmatrix_IMR90_core <- AUCmatrix_IMR90[unlist(strsplit(rownames(AUCmatrix_IMR90), "_"))[c(T,F,F,F)] %in% IMR90_HepG2_preprocess_stats_core$TFs,]



AUCmatrix_HepG2 <- matrix(0, nrow=0, ncol=12)
colnames(AUCmatrix_HepG2) <- colnames(AUCmatrix_IMR90)
cell_line <- "HEPG2"
for(DNAaccess in c("ATAC", "DNase", "MNase")){
  for(validation_regions in c("topRegions", "middleRegions", "bottomRegions")){
    load(paste0("objects_Ale/AUCs_",validation_regions,"_",cell_line,"_",DNAaccess,"_allQDAs.RData"))
    rownames(mymat_top) <- paste0(rownames(mymat_top),"_",DNAaccess,"_HepG2")
    AUCmatrix_HepG2 <- rbind(AUCmatrix_HepG2, mymat_top)
  }
}

AUCmatrix_HepG2_core <- AUCmatrix_HepG2[unlist(strsplit(rownames(AUCmatrix_HepG2), "_"))[c(T,F,F,F)] %in% IMR90_HepG2_preprocess_stats_core$TFs,]

AUCmatrix_IMR90_HepG2 <- rbind(AUCmatrix_IMR90, AUCmatrix_HepG2)

AUCmatrix_IMR90_HepG2_core <- AUCmatrix_IMR90_HepG2[unlist(strsplit(rownames(AUCmatrix_IMR90_HepG2), "_"))[c(T,F,F,F)] %in% IMR90_HepG2_preprocess_stats_core$TFs,]

AUCmatrix_IMR90_HepG2_core_clean <- AUCmatrix_IMR90_HepG2_core
rownames(AUCmatrix_IMR90_HepG2_core_clean) <- gsub("_", " ", rownames(AUCmatrix_IMR90_HepG2_core_clean))
rownames(AUCmatrix_IMR90_HepG2_core_clean) <- gsub("bottomRegions", "weak", rownames(AUCmatrix_IMR90_HepG2_core_clean))
rownames(AUCmatrix_IMR90_HepG2_core_clean) <- gsub("middleRegions", "medium", rownames(AUCmatrix_IMR90_HepG2_core_clean))
rownames(AUCmatrix_IMR90_HepG2_core_clean) <- gsub("topRegions", "strong", rownames(AUCmatrix_IMR90_HepG2_core_clean))


write.table(AUCmatrix_IMR90_HepG2_core_clean, file="data/Table_S9_core.csv", sep=",", row.names = TRUE)


pdf("Figure_S6_maxAUC_IMR90_HepG2.pdf", width=5, height=4,pointsize = 12)
m <- hist(apply(AUCmatrix_IMR90_HepG2_core,1,max ), col=cbbPalette[3],
          breaks=seq(0,1, by=0.1), main="Maximum AUC", ylab="count", xlab="AUC", ylim=c(0,80))
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)
dev.off()

sum(apply(AUCmatrix_IMR90_HepG2_core,1,max ) > 0.8)/nrow(AUCmatrix_IMR90_HepG2_core)

h_top <- hist(apply(AUCmatrix_IMR90_HepG2_core[grep("top",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ), breaks=seq(0,1, by=0.1), plot=FALSE)
h_middle <- hist(apply(AUCmatrix_IMR90_HepG2_core[grep("middle",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ), breaks=seq(0,1, by=0.1), plot=FALSE)
h_bottom <- hist(apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ), breaks=seq(0,1, by=0.1), plot=FALSE)



d_top <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("top",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))
d_middle <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("middle",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))
d_bottom <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))

wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("top",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
            apply(AUCmatrix_IMR90_HepG2_core[grep("middle",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("top",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
            apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("middle",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
            apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)


d_ATAC <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("ATAC",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))
d_DNase <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("DNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))
d_MNase <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("MNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))
d_NOMe <- density(apply(AUCmatrix_IMR90_HepG2_core[grep("NOMe",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))



wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("ATAC",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
            apply(AUCmatrix_IMR90_HepG2_core[grep("DNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("ATAC",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("MNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("ATAC",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("NOMe",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("DNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("MNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("DNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("NOMe",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("MNase",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
            apply(AUCmatrix_IMR90_HepG2_core[grep("NOMe",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("top",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)

formatC(wilcox.test(apply(AUCmatrix_IMR90_HepG2_core[grep("middle",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ),
                    apply(AUCmatrix_IMR90_HepG2_core[grep("bottom",rownames(AUCmatrix_IMR90_HepG2_core)), ],1,max ))$p.value, format = "e", digits = 2)



pdf("Figure_5_maxAUC_IMR90_HepG2_regions.pdf", width=12, height=4,pointsize = 12)
par(cex=1.2)
par(las=1)
par(mfrow=c(1,3))

m <- hist(apply(AUCmatrix_IMR90_HepG2_core,1,max ), col=cbbPalette[3],
          breaks=seq(0,1, by=0.1), main="Maximum AUC", ylab="count", xlab="AUC", ylim=c(0,70))
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)

# hist(apply(AUCmatrix_IMR90_HepG2[grep("top",rownames(AUCmatrix_IMR90_HepG2)), ],1,max ),
#      breaks=seq(0,1, by=0.1), col=rgb(0,114,178,125, maxColorValue = 255),
#      main="Maximum AUC", ylab="count", xlab="AUC", ylim=c(0,50))
# hist(apply(AUCmatrix_IMR90_HepG2[grep("middle",rownames(AUCmatrix_IMR90_HepG2)), ],1,max ),
#      breaks=seq(0,1, by=0.1), col=rgb(204,121,167,125, maxColorValue = 255), add = TRUE)
# hist(apply(AUCmatrix_IMR90_HepG2[grep("bottom",rownames(AUCmatrix_IMR90_HepG2)), ],1,max ),
#      breaks=seq(0,1, by=0.1), col=rgb(213,94,0,125, maxColorValue = 255), add = TRUE)
# mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)


plot(d_top, col=cbbPalette[6], type="l", lty=1, lwd=2,
     main="Maximum AUC (ChIP signal strength)", ylab="density", xlab="AUC", xlim=c(0,1.05))
lines(d_middle, col=cbbPalette[8], lty=1, lwd=2)
lines(d_bottom, col=cbbPalette[7], lty=1, lwd=2)
legend("topleft", col=cbbPalette[c(6,8,7)], lty=1, lwd=2,
       legend=c("strong", "medium", "weak"), bty="n", cex=1.0, horiz = F)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)


plot(d_DNase, col=cbbPalette[3], type="l", lty=1, lwd=2,
     main="Maximum AUC (accessibility method)", ylab="density", xlab="AUC", xlim=c(0,1.05), ylim=c(0,6))
lines(d_ATAC, col=cbbPalette[4], lty=1, lwd=2)
lines(d_MNase, col=cbbPalette[2], lty=1, lwd=2)
lines(d_NOMe, col=cbbPalette[5], lty=1, lwd=2)
legend("topleft", col=cbbPalette[c(3,4,2,5)], lty=1, lwd=2,
       legend=c("DNase-seq", "ATAC-seq", "MNase-seq", "NOMe-seq"), bty="n", cex=1.0, horiz = F)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)


dev.off()


load("objects_Ale/accessibilityLevels_allQDAs_allMethods_IMR90.Rda")
QDA_thresholds_IMR90 <- bigM

load("objects_Ale/accessibilityLevels_allQDAs_allMethods_HepG2.Rda")
QDA_thresholds_HepG2 <- bigM
#load("objects_Ale/accessibilityLevels_allQDAs_MNase_HepG2.Rda")
#QDA_thresholds_HepG2[,3] <- mymat3

pdf("Figure_S5_DNA_accessibility_thresholds.pdf", width=15, height=6,pointsize = 12)
par(cex=1.2)
par(las=1)
par(mfrow=c(1,2))


plot(as.numeric(row.names(QDA_thresholds_IMR90)), 100*QDA_thresholds_IMR90[,2],
     col=cbbPalette[3], type="l", lty=1, lwd=2,
     main="IMR90", ylab="% of accessible DNA", xlab="QDA")
lines(as.numeric(row.names(QDA_thresholds_IMR90)), 100*QDA_thresholds_IMR90[,1],
      col=cbbPalette[4], lty=1, lwd=2)
lines(as.numeric(row.names(QDA_thresholds_IMR90)), 100*QDA_thresholds_IMR90[,3],
      col=cbbPalette[2], lty=1, lwd=2)
lines(as.numeric(row.names(QDA_thresholds_IMR90)), 100*QDA_thresholds_IMR90[,4],
      col=cbbPalette[5], lty=1, lwd=2)
legend("topright", col=cbbPalette[c(3,4,2,5)], lty=1, lwd=2,
       legend=c("DNase-seq", "ATAC-seq", "MNase-seq", "NOMe-seq"), bty="n", cex=1.0, horiz = F)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)


plot(as.numeric(row.names(QDA_thresholds_HepG2)), 100*QDA_thresholds_HepG2[,2],
     col=cbbPalette[3], type="l", lty=1, lwd=2,
     main="HepG2", ylab="% of accessible DNA", xlab="QDA")
lines(as.numeric(row.names(QDA_thresholds_HepG2)), 100*QDA_thresholds_HepG2[,1],
      col=cbbPalette[4], lty=1, lwd=2)
lines(as.numeric(row.names(QDA_thresholds_HepG2)), 100*QDA_thresholds_HepG2[,3],
      col=cbbPalette[2], lty=1, lwd=2)

legend("topright", col=cbbPalette[c(3,4,2)], lty=1, lwd=2,
       legend=c("DNase-seq", "ATAC-seq", "MNase-seq"), bty="n", cex=1.0, horiz = F)
mtext(LETTERS[4], line = 0.7, adj = -0.09, cex=1.6)


dev.off()


quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
QDA_ID <- 1:12

#
TF_high <- apply(AUCmatrix_IMR90_HepG2_core[,9:11],1,mean, na.rm=TRUE)
TF_low <- apply(AUCmatrix_IMR90_HepG2_core[,1:3],1,mean, na.rm=TRUE)
TF_diff <- TF_high - TF_low
manual_class <- rep("other", nrow(AUCmatrix_IMR90_HepG2_core))
manual_class[which((TF_low >= 0.8 | TF_high >= 0.8) & abs(TF_diff) < 0.1)] <- "AIF"
manual_class[which((TF_high >= 0.8) & TF_diff >= 0.3)] <- "ADF"
manual_class[which((TF_low >= 0.8) & TF_diff <= -0.3)] <- "IDF"
manual_class[which(!(TF_low >= 0.65 | TF_high >= 0.65))] <- "poorly predicted"
manual_class[which((TF_high >= 0.65) & (TF_diff >= 0.1 & TF_diff < 0.3))] <- "partial AIF/ADF"
manual_class[which((TF_low >= 0.65) & (TF_diff < -0.1 & TF_diff >= -0.3))] <- "partial AIF/IDF"
table(manual_class)

AUCmatrix_IMR90_HepG2_core_manual_class <- as.data.frame(AUCmatrix_IMR90_HepG2_core)
rownames(AUCmatrix_IMR90_HepG2_core_manual_class) <- gsub("_", " ", rownames(AUCmatrix_IMR90_HepG2_core_manual_class))
rownames(AUCmatrix_IMR90_HepG2_core_manual_class) <- gsub("bottomRegions", "weak", rownames(AUCmatrix_IMR90_HepG2_core_manual_class))
rownames(AUCmatrix_IMR90_HepG2_core_manual_class) <- gsub("middleRegions", "medium", rownames(AUCmatrix_IMR90_HepG2_core_manual_class))
rownames(AUCmatrix_IMR90_HepG2_core_manual_class) <- gsub("topRegions", "strong", rownames(AUCmatrix_IMR90_HepG2_core_manual_class))
#rownames(AUCmatrix_IMR90_HepG2_core_manual_class) <- gsub("IMR90", "", rownames(AUCmatrix_IMR90_HepG2_core_manual_class))
AUCmatrix_IMR90_HepG2_core_manual_class <- split(AUCmatrix_IMR90_HepG2_core_manual_class, manual_class)


colfunc<-colorRampPalette(c("white",cbbPalette[6]))
Colors <- colfunc(20)

# all
pos <- 1:12
pdf(file = paste0("Kmeans_IMR90_threshold_heatmap_large.pdf"), width=9, height=20,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:2, ncol=1),
       height=c(nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["AIF"]])/nrow(AUCmatrix_IMR90),
                nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["poorly predicted"]])/nrow(AUCmatrix_IMR90)))
par(mar=c(4, 17.0, 3.5, 0)+0.1);

# cluster 3
#clust <- AUCmatrix_IMR90_manual_class[["AIF"]][order(rownames(AUCmatrix_IMR90_manual_class[["AIF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_manual_class[["AIF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_manual_class[["AIF"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["AIF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("AIFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[3], line = 0.7, adj = -0.09, cex=1.6)



# cluster 6
#clust <- AUCmatrix_IMR90_manual_class[["poorly predicted"]][order(rownames(AUCmatrix_IMR90_manual_class[["poorly predicted"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_manual_class[["poorly predicted"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_manual_class[["poorly predicted"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["poorly predicted"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("poorly predicted (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[6], line = 0.7, adj = -0.09, cex=1.6)


dev.off()


pdf(file = paste0("Kmeans_IMR90_threshold_heatmap_small.pdf"), width=9, height=20,pointsize = 12)
par(cex=1.1);
par(las=1)
options(scipen=999)

layout(matrix(1:4, ncol=1),
       height=c(nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]])*0.65/nrow(AUCmatrix_IMR90),
                nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]])/nrow(AUCmatrix_IMR90),
                nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]])*0.65/nrow(AUCmatrix_IMR90),
                nrow(AUCmatrix_IMR90_HepG2_core_manual_class[["other"]])*0.65/nrow(AUCmatrix_IMR90)))
par(mar=c(4, 17.0, 3.5, 0)+0.1);


# cluster 2
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]][order(rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[2], line = 0.7, adj = -0.09, cex=1.6)


# cluster 1
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]][order(rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("ADFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[1], line = 0.7, adj = -0.09, cex=1.6)



# cluster 4
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]][order(rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("partial AIFs/IDFs (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[4], line = 0.7, adj = -0.09, cex=1.6)

# cluster 5
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["other"]][order(rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["other"]])),pos]
#data.dendro <- as.dendrogram(hclust(dist(x = AUCmatrix_IMR90_HepG2_core_manual_class[["other"]][,pos])))
#data.order <- order.dendrogram(data.dendro)
#clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["other"]][data.order,pos]
clust <- AUCmatrix_IMR90_HepG2_core_manual_class[["other"]][,pos]

image(1:ncol(clust),1:nrow(clust),t(clust),
      axes = FALSE, xlab=" ", ylab=" ",col=Colors, zlim=c(0.5,1))
title(main=paste0("others (",nrow(clust),")"),cex.main=1.8)
#title( ylab="TFs", line=3.5, cex.lab=1.2, las=2)
title(xlab="QDAs",line=4.5, cex.lab=1.2)
axis(1, at=1:ncol(clust), labels=quantVec[pos],las = 1,cex.axis=1.2)
axis(LEFT <-2, at=1:nrow(clust), labels=rownames(clust),las= HORIZONTAL<-1,cex.axis=1.2)
mtext(LETTERS[5], line = 0.7, adj = -0.09, cex=1.6)

dev.off()


AUCmatrix_IMR90_HepG2_core_annotation <- matrix(unlist(strsplit(rownames(AUCmatrix_IMR90_HepG2_core), "_")), ncol=4, byrow = TRUE)
AUCmatrix_IMR90_HepG2_core_annotation <- cbind(AUCmatrix_IMR90_HepG2_core_annotation, manual_class)
AUCmatrix_IMR90_HepG2_core_annotation[,2] <- gsub("topRegions", "strong", AUCmatrix_IMR90_HepG2_core_annotation[,2])
AUCmatrix_IMR90_HepG2_core_annotation[,2] <- gsub("middleRegions", "medium", AUCmatrix_IMR90_HepG2_core_annotation[,2])
AUCmatrix_IMR90_HepG2_core_annotation[,2] <- gsub("bottomRegions", "weak", AUCmatrix_IMR90_HepG2_core_annotation[,2])

cols <- c(cbPalette[c(6,4,5,8,2,1)])

class_levels <- rownames(mm10_K562_counts_manual)
class_levels <- gsub("AIFs","AIF", class_levels)
class_levels <- gsub("ADFs","ADF", class_levels)
class_levels <- gsub("IDFs","IDF", class_levels)
class_levels <- gsub("others","other", class_levels)


AUCmatrix_IMR90_HepG2_core_annotation_df <- as.data.frame(AUCmatrix_IMR90_HepG2_core_annotation)
AUCmatrix_IMR90_HepG2_core_annotation_df[,5] <- factor(AUCmatrix_IMR90_HepG2_core_annotation[,5], levels = class_levels)
colnames(AUCmatrix_IMR90_HepG2_core_annotation_df) <- c("TF", "region", "method", "cell", "class")
AUCmatrix_IMR90_HepG2_core_annotation_df$region <- factor(AUCmatrix_IMR90_HepG2_core_annotation_df$region, levels=c("strong", "medium", "weak"))
AUCmatrix_IMR90_HepG2_core_annotation_df <- AUCmatrix_IMR90_HepG2_core_annotation_df[order(AUCmatrix_IMR90_HepG2_core_annotation_df$method),]
AUCmatrix_IMR90_HepG2_core_annotation_df <- AUCmatrix_IMR90_HepG2_core_annotation_df[order(AUCmatrix_IMR90_HepG2_core_annotation_df$region),]
AUCmatrix_IMR90_HepG2_core_annotation_df <- AUCmatrix_IMR90_HepG2_core_annotation_df[order(AUCmatrix_IMR90_HepG2_core_annotation_df$TF),]


AUCmatrix_IMR90_HepG2_core_annotation_df_cell <- split(AUCmatrix_IMR90_HepG2_core_annotation_df, AUCmatrix_IMR90_HepG2_core_annotation_df$cell)

for(i in 1:length(AUCmatrix_IMR90_HepG2_core_annotation_df_cell)){
  regions <- split(AUCmatrix_IMR90_HepG2_core_annotation_df_cell[[i]], AUCmatrix_IMR90_HepG2_core_annotation_df_cell[[i]]$region)

  buffer_merged <- NULL
  for(j in 1:length(regions)){
    buffer <- as.data.frame(matrix(regions[[j]][,5], ncol=length(unique(regions[[j]][,3])), byrow = TRUE))
    for(k in 1:ncol(buffer)){
      buffer[,k] <- as.numeric(factor(buffer[,k], levels=class_levels))
    }
    colnames(buffer) <- unique(regions[[j]][,3])
    rownames(buffer) <- unique(regions[[j]][,1])
    buffer_extra <- matrix(NA, ncol=length(unique(regions[[j]][,3])), nrow=(11-nrow(buffer)))
    if(names(AUCmatrix_IMR90_HepG2_core_annotation_df_cell)[i] == "IMR90" & nrow(buffer_extra) == 1){
      rownames(buffer_extra) <- "ELK1"
    }
    colnames(buffer_extra) <- colnames(buffer)
    buffer <- rbind(buffer, buffer_extra)
    if(ncol(buffer) == 3){
      buffer <- cbind(buffer, NA)
    }
    buffer <- buffer[order(rownames(buffer), decreasing = T),]

    if(is.null(buffer_merged)){
      buffer_merged <- buffer
    } else{
      buffer_merged <- cbind(buffer_merged, buffer)
    }
  }
#  ids<-which(buffer_merged > 4,arr.ind = TRUE)
#  for(l in 1:nrow(ids)){
#    buffer_merged[ids[l,1],ids[l,2]] <- NA
#  }
  order_id <- order(apply(buffer_merged,1,mean, na.rm=TRUE), decreasing = T)

  for(j in 1:length(regions)){
    buffer <- as.data.frame(matrix(regions[[j]][,5], ncol=length(unique(regions[[j]][,3])), byrow = TRUE))
    for(k in 1:ncol(buffer)){
      buffer[,k] <- as.numeric(factor(buffer[,k], levels=class_levels))
    }
    colnames(buffer) <- unique(regions[[j]][,3])
    rownames(buffer) <- unique(regions[[j]][,1])
    buffer_extra <- matrix(NA, ncol=length(unique(regions[[j]][,3])), nrow=(11-nrow(buffer)))
    if(names(AUCmatrix_IMR90_HepG2_core_annotation_df_cell)[i] == "IMR90" & nrow(buffer_extra) == 1){
      rownames(buffer_extra) <- "ELK1"
    }
    colnames(buffer_extra) <- colnames(buffer)
    buffer <- rbind(buffer, buffer_extra)
    if(ncol(buffer) == 3){
      buffer <- cbind(buffer, NA)
    }
    buffer <- buffer[order(rownames(buffer), decreasing = T),]

    buffer <- buffer[order_id,]

    pdf(file = paste0("summary_analysis_",names(AUCmatrix_IMR90_HepG2_core_annotation_df_cell)[i],"_",names(regions)[j],"_heatmap.pdf"),
        width=7, height=(2.65+nrow(buffer)*0.3),pointsize = 12)
    par(cex=1.1);
    par(las=1)
    options(scipen=999)

    #layout(matrix(1:2, ncol=2), width=c(1,0.7))
    par(mar=c(5, 9.0, 3.5, 13)+0.1, xpd=TRUE);
    image(1:ncol(buffer),1:nrow(buffer),t(buffer), axes = FALSE,
          xlab=" ", ylab=" ",col=cols, zlim=c(1,length(class_levels)), main=paste0(names(regions)[j], " regions"))
    segments(0.5,0:nrow(buffer)+0.5, ncol(buffer)+0.5,0:nrow(buffer)+0.5, lty=1, lwd=2)
    segments(0:ncol(buffer)+0.5, 0.5,0:ncol(buffer)+0.5,nrow(buffer)+0.5, lty=1, lwd=2)
    axis(1, at=1:ncol(buffer), labels=colnames(buffer),las = 2,cex.axis=1.2, tick = F)
    axis(LEFT <-2, at=1:nrow(buffer), labels=rownames(buffer),las= HORIZONTAL<-1,cex.axis=1.2, tick = F)
    #plot(0,0,type="n")
    legend("right", inset=c(-1,0),fill=cols, legend = class_levels, bty="n")
    dev.off()
  }
}





IMR90_HepG2_counts_manual <- matrix(0, ncol=length(unique(AUCmatrix_IMR90_HepG2_core_annotation[,1])), nrow=6)
rownames(IMR90_HepG2_counts_manual) <- c("AIFs", "partial AIFs/ADFs", "ADFs", "partial AIFs/IDFs", "others", "poorly predicted")
colnames(IMR90_HepG2_counts_manual) <- unique(AUCmatrix_IMR90_HepG2_core_annotation[,1])

for(i in 1:ncol(IMR90_HepG2_counts_manual)){
  TF <- colnames(IMR90_HepG2_counts_manual)[i]
  IMR90_HepG2_counts_manual[1,i] <- IMR90_HepG2_counts_manual[1,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["AIF"]])))
  IMR90_HepG2_counts_manual[2,i] <- IMR90_HepG2_counts_manual[2,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/ADF"]])))
  IMR90_HepG2_counts_manual[3,i] <- IMR90_HepG2_counts_manual[3,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["ADF"]])))
  IMR90_HepG2_counts_manual[4,i] <- IMR90_HepG2_counts_manual[4,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["partial AIF/IDF"]])))
  IMR90_HepG2_counts_manual[5,i] <- IMR90_HepG2_counts_manual[5,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["other"]])))
  IMR90_HepG2_counts_manual[6,i] <- IMR90_HepG2_counts_manual[6,i] +
    length(grep(paste0(TF," "), rownames(AUCmatrix_IMR90_HepG2_core_manual_class[["poorly predicted"]])))
}

pdf("Figure_5_classifications.pdf", width=8, height=4,pointsize = 12)
cols <- c(cbPalette[c(6,4,3,8,1,2)])
#cols <- c(cbPalette[c(7,5,4,2,3,1)])
cols <- c(cbPalette[c(6,5,4,8,2,1)])
cols <- c(cbPalette[c(6,4,5,8,2,1)])
par(las=2)
barplot(IMR90_HepG2_counts_manual, col=cols, main="classifications of TFs in IMR90 and HepG2", ylab="count", xlim=c(0, 25))
legend("topright", fill=cols, legend = rownames(IMR90_HepG2_counts_manual), bty="n", horiz = F)
dev.off()


################################################################################
# TF behavior in Her2 over expression in MCF10A
################################################################################

################################################################################
# JUN
################################################################################

JUN_profiles <- list("")
JUN_profiles_mean <- list("")
load("objects_Ale/Her2OE/JUN/new/predicted_chip_ctrl_newRegions_MCF10_JUN_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
JUN_profiles[["gained_ctrl"]] <- buffer
JUN_profiles_mean[["gained_ctrl"]] <- apply(JUN_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/JUN/new/predicted_chip_her2_newRegions_MCF10_JUN_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
JUN_profiles[["gained_her2oe"]] <- buffer
JUN_profiles_mean[["gained_her2oe"]] <- apply(JUN_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/JUN/lost/predicted_chip_ctrl_lostRegions_MCF10_JUN_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
JUN_profiles[["lost_ctrl"]] <- buffer
JUN_profiles_mean[["lost_ctrl"]] <- apply(JUN_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/JUN/lost/predicted_chip_her2_lostRegions_MCF10_JUN_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
JUN_profiles[["lost_her2oe"]] <- buffer
JUN_profiles_mean[["lost_her2oe"]] <- apply(JUN_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(JUN_profiles_mean, JUN_profiles, file="objects_Ale/Her2OE/JUN_profiles_mean.RData")

################################################################################
# JUND
################################################################################

JUND_profiles <- list("")
JUND_profiles_mean <- list("")
load("objects_Ale/Her2OE/JUND/new/predicted_chip_ctrl_newRegions_MCF10_JUND_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
JUND_profiles[["gained_ctrl"]] <- buffer
JUND_profiles_mean[["gained_ctrl"]] <- apply(JUND_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/JUND/new/predicted_chip_her2_newRegions_MCF10_JUND_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
JUND_profiles[["gained_her2oe"]] <- buffer
JUND_profiles_mean[["gained_her2oe"]] <- apply(JUND_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/JUND/lost/predicted_chip_ctrl_lostRegions_MCF10_JUND_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
JUND_profiles[["lost_ctrl"]] <- buffer
JUND_profiles_mean[["lost_ctrl"]] <- apply(JUND_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/JUND/lost/predicted_chip_her2_lostRegions_MCF10_JUND_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
JUND_profiles[["lost_her2oe"]] <- buffer
JUND_profiles_mean[["lost_her2oe"]] <- apply(JUND_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(JUND_profiles_mean, JUND_profiles, file="objects_Ale/Her2OE/JUND_profiles_mean.RData")

################################################################################
# ATF1
################################################################################

ATF1_profiles <- list("")
ATF1_profiles_mean <- list("")
load("objects_Ale/Her2OE/ATF1/new/predicted_chip_ctrl_newRegions_MCF10_ATF1_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
ATF1_profiles[["gained_ctrl"]] <- buffer
ATF1_profiles_mean[["gained_ctrl"]] <- apply(ATF1_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/ATF1/new/predicted_chip_her2_newRegions_MCF10_ATF1_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
ATF1_profiles[["gained_her2oe"]] <- buffer
ATF1_profiles_mean[["gained_her2oe"]] <- apply(ATF1_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/ATF1/lost/predicted_chip_ctrl_lostRegions_MCF10_ATF1_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
ATF1_profiles[["lost_ctrl"]] <- buffer
ATF1_profiles_mean[["lost_ctrl"]] <- apply(ATF1_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/ATF1/lost/predicted_chip_her2_lostRegions_MCF10_ATF1_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
ATF1_profiles[["lost_her2oe"]] <- buffer
ATF1_profiles_mean[["lost_her2oe"]] <- apply(ATF1_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(ATF1_profiles_mean, ATF1_profiles, file="objects_Ale/Her2OE/ATF1_profiles_mean.RData")

################################################################################
# ETV6
################################################################################

ETV6_profiles <- list("")
ETV6_profiles_mean <- list("")
load("objects_Ale/Her2OE/ETV6/new/predicted_chip_ctrl_newRegions_MCF10_ETV6_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
ETV6_profiles[["gained_ctrl"]] <- buffer
ETV6_profiles_mean[["gained_ctrl"]] <- apply(ETV6_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/ETV6/new/predicted_chip_her2_newRegions_MCF10_ETV6_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
ETV6_profiles[["gained_her2oe"]] <- buffer
ETV6_profiles_mean[["gained_her2oe"]] <- apply(ETV6_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/ETV6/lost/predicted_chip_ctrl_lostRegions_MCF10_ETV6_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
ETV6_profiles[["lost_ctrl"]] <- buffer
ETV6_profiles_mean[["lost_ctrl"]] <- apply(ETV6_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/ETV6/lost/predicted_chip_her2_lostRegions_MCF10_ETV6_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
ETV6_profiles[["lost_her2oe"]] <- buffer
ETV6_profiles_mean[["lost_her2oe"]] <- apply(ETV6_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(ETV6_profiles_mean, ETV6_profiles, file="objects_Ale/Her2OE/ETV6_profiles_mean.RData")

################################################################################
# MYC
################################################################################

MYC_profiles <- list("")
MYC_profiles_mean <- list("")
load("objects_Ale/Her2OE/MYC/new/predicted_chip_ctrl_newRegions_MCF10_MYC_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
MYC_profiles[["gained_ctrl"]] <- buffer
MYC_profiles_mean[["gained_ctrl"]] <- apply(MYC_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/MYC/new/predicted_chip_her2_newRegions_MCF10_MYC_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
MYC_profiles[["gained_her2oe"]] <- buffer
MYC_profiles_mean[["gained_her2oe"]] <- apply(MYC_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/MYC/lost/predicted_chip_ctrl_lostRegions_MCF10_MYC_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
MYC_profiles[["lost_ctrl"]] <- buffer
MYC_profiles_mean[["lost_ctrl"]] <- apply(MYC_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/MYC/lost/predicted_chip_her2_lostRegions_MCF10_MYC_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
MYC_profiles[["lost_her2oe"]] <- buffer
MYC_profiles_mean[["lost_her2oe"]] <- apply(MYC_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(MYC_profiles_mean, MYC_profiles, file="objects_Ale/Her2OE/MYC_profiles_mean.RData")

################################################################################
# NFATC3
################################################################################

NFATC3_profiles <- list("")
NFATC3_profiles_mean <- list("")
load("objects_Ale/Her2OE/NFATC3/new/predicted_chip_ctrl_newRegions_MCF10_NFATC3_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
NFATC3_profiles[["gained_ctrl"]] <- buffer
NFATC3_profiles_mean[["gained_ctrl"]] <- apply(NFATC3_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/NFATC3/new/predicted_chip_her2_newRegions_MCF10_NFATC3_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
NFATC3_profiles[["gained_her2oe"]] <- buffer
NFATC3_profiles_mean[["gained_her2oe"]] <- apply(NFATC3_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/NFATC3/lost/predicted_chip_ctrl_lostRegions_MCF10_NFATC3_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
NFATC3_profiles[["lost_ctrl"]] <- buffer
NFATC3_profiles_mean[["lost_ctrl"]] <- apply(NFATC3_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/NFATC3/lost/predicted_chip_her2_lostRegions_MCF10_NFATC3_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
NFATC3_profiles[["lost_her2oe"]] <- buffer
NFATC3_profiles_mean[["lost_her2oe"]] <- apply(NFATC3_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(NFATC3_profiles_mean, NFATC3_profiles, file="objects_Ale/Her2OE/NFATC3_profiles_mean.RData")

################################################################################
# SRF
################################################################################

SRF_profiles <- list("")
SRF_profiles_mean <- list("")
load("objects_Ale/Her2OE/SRF/new/predicted_chip_ctrl_newRegions_MCF10_SRF_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
SRF_profiles[["gained_ctrl"]] <- buffer
SRF_profiles_mean[["gained_ctrl"]] <- apply(SRF_profiles[["gained_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/SRF/new/predicted_chip_her2_newRegions_MCF10_SRF_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
SRF_profiles[["gained_her2oe"]] <- buffer
SRF_profiles_mean[["gained_her2oe"]] <- apply(SRF_profiles[["gained_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)

load("objects_Ale/Her2OE/SRF/lost/predicted_chip_ctrl_lostRegions_MCF10_SRF_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_ctrl@profiles[[1]])){
  print(paste0(i,"/",length(chip_ctrl@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_ctrl@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_ctrl@profiles[[1]][[i]]$ChIP)
  }
}
SRF_profiles[["lost_ctrl"]] <- buffer
SRF_profiles_mean[["lost_ctrl"]] <- apply(SRF_profiles[["lost_ctrl"]],2,mean)

rm(buffer)
rm(chip_ctrl)

load("objects_Ale/Her2OE/SRF/lost/predicted_chip_her2_lostRegions_MCF10_SRF_gw.Rda")
buffer <- NULL
for(i in 1:length(chip_her2@profiles[[1]])){
  print(paste0(i,"/",length(chip_her2@profiles[[1]])))
  if(is.null(buffer)){
    buffer <- rbind(chip_her2@profiles[[1]][[i]]$ChIP)
  } else{
    buffer <- rbind(buffer, chip_her2@profiles[[1]][[i]]$ChIP)
  }
}
SRF_profiles[["lost_her2oe"]] <- buffer
SRF_profiles_mean[["lost_her2oe"]] <- apply(SRF_profiles[["lost_her2oe"]],2,mean)
rm(buffer)
rm(chip_her2)
save(SRF_profiles_mean, SRF_profiles, file="objects_Ale/Her2OE/SRF_profiles_mean.RData")



pdf("Figure_S7.pdf", width=25, height=4,pointsize = 12)
cols <- c(cbPalette[c(6,3,7,2)])
par(las=2)
par(mfrow=c(1,7))

buffer <- ATF1_profiles_mean
buffer[[1]] <- NULL
plot(ATF1_profiles_mean[["lost_ctrl"]], col=cols[1], main="ATF1", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(ATF1_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(ATF1_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(ATF1_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)

buffer <- ETV6_profiles_mean
buffer[[1]] <- NULL
plot(ETV6_profiles_mean[["lost_ctrl"]], col=cols[1], main="ETV6", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(ETV6_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(ETV6_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(ETV6_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)

buffer <- JUN_profiles_mean
buffer[[1]] <- NULL
plot(JUN_profiles_mean[["lost_ctrl"]], col=cols[1], main="JUN", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(JUN_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(JUN_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(JUN_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)

buffer <- JUND_profiles_mean
buffer[[1]] <- NULL
plot(JUND_profiles_mean[["lost_ctrl"]], col=cols[1], main="JUND", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(JUND_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(JUND_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(JUND_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)


buffer <- MYC_profiles_mean
buffer[[1]] <- NULL
plot(MYC_profiles_mean[["lost_ctrl"]], col=cols[1], main="MYC", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(MYC_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(MYC_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(MYC_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)

buffer <- NFATC3_profiles_mean
buffer[[1]] <- NULL
plot(NFATC3_profiles_mean[["lost_ctrl"]], col=cols[1], main="NFATC3", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(NFATC3_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(NFATC3_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(NFATC3_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)


buffer <- SRF_profiles_mean
buffer[[1]] <- NULL
plot(SRF_profiles_mean[["lost_ctrl"]], col=cols[1], main="SRF", ylab="", type="l", lty=1, lwd=2,
     ylim=c(min(unlist(buffer)), 1.3*max(unlist(buffer))), xaxt="n", xlab="")
lines(SRF_profiles_mean[["lost_her2oe"]], col=cols[2], lty=1, lwd=2)
lines(SRF_profiles_mean[["gained_ctrl"]], col=cols[3], lty=1, lwd=2)
lines(SRF_profiles_mean[["gained_her2oe"]], col=cols[4], lty=1, lwd=2)
axis(BELOW<-1, at=c(0,22,44), labels=c("-20", "0", "20"), las=1)
legend("topright", col=cols, lty=1, legend = c("lost (control)", "lost (Her2)","gained (control)", "gained (Her2)"), lwd=2, bty="n", horiz = F, cex = 1.0)


dev.off()
