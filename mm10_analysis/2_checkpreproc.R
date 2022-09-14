
################################################################################
# libraries and functions
################################################################################
setwd("/storage/projects/ZabetLab/mm10_ENCODE/data/")

library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(MotifDb)
library(ChIPanalyser)

peaksLoading<-function(x){

    if(grepl(x=x,pattern=".bed")& !grepl(x=x,pattern=".gff")){

          x<-read.table(x, stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }

    }else if(grepl(x=x, pattern="_gw.Rda")){
        x<-get(load(x))
    }else if(grepl(x=x, pattern=".gff3")){
          x<-read.table(x, skip=30,stringsAsFactors=F)
          if(length(grep(x=x[,1], pattern="chr"))==0){
              x[,1]<-paste0("chr",x[,1])
              x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
          }else{
              x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
          }

    } else if(grepl(x=x, pattern=".gff")){
        x<-read.table(x, skip=30,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }
    } else{
        x<-read.table(x, header=T,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
             x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }
    }
    return(x)
}

#check the number of peaks and sum width of peaks for any preprocessed TF
howManyPeaks<- function (TF){
	peaks_data <- peaksLoading(paste0("../data/",TF,"_CH12_peaks_peaks.narrowPeak"))
	x<- length(peaks_data)
	y<-sum(peaks_data@ranges@width)
	result<- cat(TF, "Number of peaks:", x,",total width of peaks:", y,"\n", fill=TRUE)
	return (result)
}

#get the number of peaks
peaksPlot<- function(TF){
	peaks_data <- peaksLoading(paste0("../data/",TF,"_",cell_line,"_peaks_peaks.narrowPeak"))
	x<- length(peaks_data)
}

###################
##get values
###################
#identify the cell line
cell_line<-"CH12"

#apply the peaksPlot function to a list of TFs
#make a table containing TF names and number of ChIPseq peaks
names = readLines(paste0("../scripts/TFs_list_", cell_line,".txt"))
peaks = sapply(names, peaksPlot)
df<-rbind(names,peaks)
write.table(peaks, paste0("../scripts/numberOfPeaks_gw_TFs_",cell_line,".txt"))

#plot number of ChIPseq peaks
pdf(paste0("../imgs/numberOfPeaks_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(peaks~names, data=df, main="Number of ChIPseq peaks",
   xlab="TFs", ylab="Number of peaks", cex.axis=0.8,yaxt="n", las=2)
      axis(2, c(0, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000), cex.axis=0.8, las=1)
dev.off()


#######
#plots
#######

#you will need a table containing the data from preprocessing (fastqc, bowtie2 and macs2)
#did this part manually
TF_table<- read.table (paste0("../data/TFs_preprocessing_data_",cell_line,".csv"), header=T , sep=",")

#number of reads from fastqc
pdf(paste0("../imgs/numberOfReads_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(reads~names, data=TF_table, main="Number of reads before trimming",
   xlab="TFs", ylab="Number of reads", cex.axis=0.8, yaxt="n",las=2)
   axis(2, c(0, 1.0e+07,2.0e+07,3.0e+07,4.0e+07,5.0e+07,6.0e+07, 7.0e+07,8.0e+07), cex.axis=0.8, las=1)
dev.off()
#number of reads from fastqc after trimming with trimmomatic
pdf(paste0("../imgs/numberOfTrimmedReads_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(trimmedreads~names, data=TF_table, main="Number of reads after trimming",
   xlab="TFs", ylab="Number of reads", cex.axis=0.8, yaxt="n",las=2)
      axis(2, c(0, 1.0e+07,2.0e+07,3.0e+07,4.0e+07,5.0e+07,6.0e+07, 7.0e+07,8.0e+07), cex.axis=0.8, las=1)
dev.off()
#number of aligned reads from bowtie2
pdf(paste0("../imgs/numberOfAlignedReads_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(aligned~names, data=TF_table, main="Number of aligned reads",
   xlab="TFs", ylab="Number of aligned reads", cex.axis=0.8,yaxt="n", las=2)
      axis(2, c(0, 1.0e+07,2.0e+07,3.0e+07,4.0e+07,5.0e+07,6.0e+07, 7.0e+07,8.0e+07), cex.axis=0.8, las=1)
dev.off()
#percent of aligned reads from bowtie2
pdf(paste0("../imgs/AlignedReadsPercent_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(alignedpercent~names, data=TF_table, ylim=c(0, 1), main="Percentage of aligned reads",
   xlab="TFs", ylab="Aligned reads (%)", cex.axis=0.8,yaxt="n", las=2)
      axis(2, c( 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1))
dev.off()
#number of peaks from macs2 (from this table - same result as before)
pdf(paste0("../imgs/numberOfPeaks_gw_barplot_TFs_",cell_line,".pdf"), width = 12, height = 10, pointsize = 12)
barplot(peaks~names, data=TF_table, main="Number of ChIPseq peaks",
   xlab="TFs", ylab="Number of peaks", cex.axis=0.8,yaxt="n", las=2)
      axis(2, c(0, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000), cex.axis=0.8, las=1)
dev.off()

#check if it has pwm


