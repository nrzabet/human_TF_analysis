################################################################
###################### Barplots - general ######################
################################################################


## Creating barplots of % alignment after Bowtie.
## screen, qrsh, cd /, R3.5.0
cell_line <- "IMR90"
setwd(paste0("storage/st20d/HiC_Kc167/alessandra/",cell_line,"/imgs/barplots/"))

## I got this info by checking the log files created after bowtie (pre processing). 
align <- sort(c(94.9, 92.55, 92.89, 94.65, 94.45, 90.35, 93.92, 96.11, 91.51, 97.09,
                92.68, 93.73, 94.57, 93.65, 89.26, 97.01, 91.60, 97.91), decreasing = TRUE)

## Make sure alignmment percentage and TF names correspond. 
tfs <- c("Ctrl 1", "CEBPB", "CHD1", "MAFK", "RCOR1", "FOS", "RAD21", "BHLHE40", "USF2", "CTCF", "POLR2A", "SMC3",
         "MAZ", "ELK1", "Ctrl 2", "RFX5", "NFE2L2", "MXI1")

pdf(paste0("alignmentPercentage.pdf"), width=20, height=16, pointsize = 13)
barplot(align,
        names.arg = tfs,
        yaxt="n", # makes y axis disappear so I can do it myself (below)
        las = 2, # makes TF names vertical
        ylim = c(0,100), # changes range of y axis
        main = paste0("Alignment % of trimmed data to reference genome after Bowtie (", cell_line," CL)"),
        col = "blue",
        ylab = "Percentage alignment")
axis(LEFT<-2,at=seq(0,100,by=20),labels=paste0(seq(0,100,by=20),"%"),las=2,cex.axis=1)
dev.off()





## Checking number of ChIP-Seq peaks.
## "Number of ChIP-Seq peaks" shows the result of MACS2 peak calling.

## screen, qrsh, R3.5.0
setwd(paste0("storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/"))

library(GenomicRanges)
## If package is not installed, type:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(GenomicRanges)

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

tfs <- c("CEBPB", "MAFK", "RAD21", "BHLHE40", "POLR2A", "FOS", "SMC3",
         "MAZ", "CTCF", "RCOR1", "MXI1", "USF2", "NFE2L2", "CHD1", "RFX5", "ELK1")

peaksPlot<- function(tfs){
        peaks_data <- peaksLoading(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/",
                                          tfs,"-human_",cell_line,"_SE_trim-peaks_peaks.narrowPeak"))
        x<- length(peaks_data)
}

peaks <- sapply(tfs, peaksPlot)
df <- rbind(tfs,peaks)
write.table(peaks, paste0("numberOfPeaks_",cell_line,"_TFs.txt"))

## Now open the file just created (numberOfPeaks_cellline_TFs.txt), 
## copy-paste the number of peaks corresponding to each TF below and put it in a variable (which here below i called chippeaks).
chippeaks <- sort(c(227924,111711,87317,65406,64179,58914,62082,57739,56686,55700,
                    46379,40012,34719,27160,5461,145), decreasing = TRUE)

## PS: make sure that the TF name and number of peaks correspond!!!!!!

setwd(paste0("storage/st20d/HiC_Kc167/alessandra/",cell_line,"/imgs/barplots/"))

pdf(paste0("numberChIP-seqPeaks.pdf"), width=20, height=16, pointsize = 13)
mybar <- barplot(chippeaks,
        names.arg = tfs,
        log = "y",
        yaxt="n",
        las = 2, #makes TF names vertical
        ylim = c(1,1000000), # changes range of y axis based on higher number
        main = paste0("Number of ChIP-seq peaks after MACS2 peak calling (",cell_line," CL)"),
        col = "red",
        ylab = "Number of peaks",
        cex.lab = 1,
        cex.axis = 0.7,
        space = c(0,0))
at <- c(1, 10, 100, 1000, 10000, 100000, 1000000)
axis(2,
    at = at,
    las=2, cex.axis=0.7)
text(mybar, 1, labels = paste("n:", chippeaks, sep=""),cex=0.7,pos=3) 
dev.off()