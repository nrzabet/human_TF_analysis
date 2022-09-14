setwd("/storage/projects/ZabetLab/mm10_ENCODE/scripts/")

cbbPalette <- c("#F79019","#FF0000", "#0000FF","#2ADCB3", "#4EA403","#9DFC4A",    "#C52ADC", "#FF1A94")


#################
#number of peaks

#have a nice table with all cell line peaks
df<- read.csv("../data/set2/peaks_all.csv", header=T, sep=",")

df
rownames(df)<-df[,1]
df<-df[,-1]
df <- df[order(df[,1], decreasing = TRUE),]

mx <- t(as.matrix(df)) 

mx_vect <- as.vector(mx)
names(mx_vect) <- paste0(rep(colnames(mx),each = nrow(mx)),"_" ,rep(rownames(mx),ncol(mx)))
mx_vect <- mx_vect[!is.na(mx_vect)]
mx_vect <- mx_vect[order(mx_vect, decreasing = TRUE)]

cells <- unlist(strsplit(names(mx_vect), "_"))[c(F,T)]
cells_unique <- unique(unlist(strsplit(names(mx_vect), "_"))[c(F,T)])

colors <- cbbPalette
names(colors) <- cells_unique
pdf("../imgs/set2/numberOfPeaks_gw_barplot_cols.pdf", width = 15, height = 5, pointsize = 12)
par(mar=c(8, 4.0, 3, 1)+0.1)
barplot(log10(mx_vect), names.arg=gsub("_","-",names(mx_vect)),  col=colors[match(cells, cells_unique)],
        main="Number of ChIPseq peaks",xlab="",yaxt="n", ylab="Number of peaks", 
        cex.axis=0.8, cex.names=0.6, las=2)
legend("topright", legend=names(colors), fill=colors, cex=0.6, bty="n")
axis(2, 0:5,c("1", "10", "100", "1000", "10,000", "100,000"), cex.axis=0.8, las=1)
dev.off()


pdf("numberOfPeaks_gw_barplot_single.pdf", width = 15, height = 5, pointsize = 12)
par(mar=c(8, 4.0, 3, 1)+0.1)
barplot(log10(mx_vect), names.arg=gsub("_","-",names(mx_vect)),  col=colors[6],
        main="Number of ChIPseq peaks",xlab="", ylab="log10(#peaks)", 
        cex.axis=0.8, cex.names=0.6, las=2)
dev.off()


###########
#alignment
df<- read.csv("../data/set2/align_all.csv", header=T, sep=",")

rownames(df)<-df[,1]

df<-df[,-1]
df <- df[order(df[,1], decreasing = TRUE),]

mx <- t(as.matrix(df))
mx <- t(as.matrix(df)) 

mx_vect <- as.vector(mx)
names(mx_vect) <- paste0(rep(colnames(mx),each = nrow(mx)),"_" ,rep(rownames(mx),ncol(mx)))
mx_vect <- mx_vect[!is.na(mx_vect)]
mx_vect <- mx_vect[order(mx_vect, decreasing = TRUE)]

cells <- unlist(strsplit(names(mx_vect), "_"))[c(F,T)]
cells_unique <- unique(unlist(strsplit(names(mx_vect), "_"))[c(F,T)])

cbbPalette <- c("#FF0000","#4EA403", "#C52ADC", "#F79019","#0000FF","#9DFC4A",    "#FF1A94","#2ADCB3")

colors <- cbbPalette
names(colors) <- cells_unique

pdf("../imgs/set2/percentOfAlignment_gw_barplot_cols.pdf", width = 15, height = 5, pointsize = 12)
par(mar=c(8, 4.0, 3, 1)+0.1)
barplot(mx_vect, names.arg=gsub("_","-",names(mx_vect)),  col=colors[match(cells, cells_unique)],
        main="Alignment percentage",xlab="", ylab="Ratio of alignment", 
        cex.axis=0.8, cex.names=0.6, las=2,ylim=c(0,1))
legend("topright", legend=names(colors), fill=colors, cex=0.6, bty="n")
dev.off()

