######### Validating results with known pioneer TFs in HepG2 with ChIPanalyser #########
## CREB1
## https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&files.file_type=fastq&biosample_ontology.term_name=HepG2&target.label=CREB1

## PREPROCESSING ##
## Open link above and select fastq only. Hit download. Open "files.txt", run first line on Google. It will download a file called "metadata".
## Change metadata name to match cell line e.g. metadata_cell_line_TF.tsv
## Move everything to cluster. Open screen, qrsh, go to right directory
xargs -L 1 curl -O -J -L < files.txt 
gunzip *.gz # if the files are zipped

R3.5.0
cell_line <- "HEPG2"
srr1 <- "ENCFF000PGX"
srr2 <- "ENCFF000PGZ"
srr3 <- "ENCFF149OFF"
srr4 <- "ENCFF199RYE"
srr5 <- "ENCFF511CNY"
srr6 <- "ENCFF384KWM"
srr7 <- "ENCFF517FOY"
srr8 <- "ENCFF951QJT" # fastq names from files

commands1<- c("#!/bin/bash\n\n")
commands1 <- paste0("cat ", srr1,".fastq ", srr2,".fastq ", srr3,".fastq ", srr4,".fastq ", srr5,".fastq ", srr6,".fastq ",
                    srr7,".fastq ", srr8,".fastq > CREB1_",cell_line,".fastq", collapse =  "\n")
write(commands1, file="download_CREB1_files_1.sh")

commands2 <- c("#!/bin/bash\n\n")
fastqc<-paste0("fastqc CREB1_",cell_line,".fastq")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U CREB1_",cell_line,".fastq -S 'CREB1_",cell_line,"_chip.sam' > 'CREB1_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t CREB1_",cell_line,"_chip.sam -f SAM -g 2.9e9 -n CREB1_",cell_line,"_peaks -B -q 0.05 > 'CREB1_",cell_line, "_peaksLog.log'")
commands2<-paste0(fastqc, "\n", bowtie, "\n", macs2, "\n")
write(commands2, file="download_CREB1_files_2.sh")

q()
n

chmod 770 download_CREB1_files_1.sh
chmod 770 download_CREB1_files_2.sh
./download_CREB1_files_1.sh
./download_CREB1_files_2.sh


########################## PRE-PROCESSING DONE! ##########################  




## FINDING NUMBER OF CHIPSEQ PEAKS ##
R3.5.0
cell_line <- "HepG"
setwd(paste0("storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/"))
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

tfs <- "CREB1"
cell_line <- "HEPG2"

peaksPlot<- function(tfs){
  peaks_data <- peaksLoading(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/",
                                    tfs,"_",cell_line,"_peaks_peaks.narrowPeak"))
  x<- length(peaks_data)
}
peaks <- sapply(tfs, peaksPlot) # number of peaks for CREB1: 205936


########################## NUMBER OF CHIP-SEQ DONE! ##########################




## GET MOTIFS ##
## screen, qrsh, cd /, go to IMR90/preProcessing folder
ls CREB1_HEPG2_peaks_treat_pileup.bdg > HEPG2_CREB1_bdg.txt
mv /storage/st20d/HiC_Kc167/alessandra/HEPG2/preProcessing/HEPG2_CREB1_bdg.txt /storage/st20d/HiC_Kc167/alessandra/HEPG2/ChIPanalyser/motifs/
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/ChIPanalyser/motifs/
R3.5.0

library(MotifDb)
cell_line <- "HEPG2"
TFnames <- "CREB1"

for(i in TFnames){
  
  if(length(query(query(query(MotifDb, "hsapiens",
                              notStrings = c("Mmusculus", "rrattus")), i), "jaspar2018"))!= 0){
    
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "jaspar2018")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    ##get the source of the PFM
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "jaspar2018"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10"))!= 0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10")
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "HOCOMOCOv10"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02")) !=0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "jolma2013")) !=0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "jolma2013")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "cisbp_1.02"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    
  }else if(length(query(query(MotifDb, "hsapiens",
                              notStrings=c("Mmusculus", "rrattus")), i))!= 0){
    PFM <- query(query(MotifDb, "hsapiens", notStrings=c("Mmusculus", "rrattus")), i)[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(MotifDb, "hsapiens"), i))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
    save(PFM, file = c(paste0(i, ".Rda")))
    
  }else if(length(query(query(query(MotifDb, "hsapiens",
                                    notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon")) != 0){
    PFM <- query(query(query(MotifDb, "hsapiens",
                             notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon")[[1]]
    save(PFM, file = c(paste0(i, ".Rda")))
    from <- matrix(c(i, names(query(query(query(MotifDb, "hsapiens",
                                                notStrings=c("Mmusculus", "rrattus")), i), "SwissRegulon"))[[1]]),
                   ncol = 2, nrow =1)
    cat(from, "\n", file = "PFMsource.csv", append = T, sep = "\t")
  }
}

########################## MOTIFS DONE! ##########################



## EXECUTING CHIPANALYSER ##
## Move ChIPanalyser_general_forKnownPioneers.R and make executable.

# ATAC, DNase, MNase and NOMe: top, middle and bottom - to do

# DNase top middle and bottom - to do
qsub -cwd -j y -q all.q -o ./t_CREB1_12.txt -b y -N t_CREB1_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 12 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_11.txt -b y -N t_CREB1_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 11 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_10.txt -b y -N t_CREB1_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 10 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_9.txt -b y -N t_CREB1_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 9 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_8.txt -b y -N t_CREB1_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 8 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_7.txt -b y -N t_CREB1_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 7 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_6.txt -b y -N t_CREB1_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 6 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_5.txt -b y -N t_CREB1_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 5 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_4.txt -b y -N t_CREB1_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 4 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_3.txt -b y -N t_CREB1_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 3 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_2.txt -b y -N t_CREB1_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 2 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_CREB1_1.txt -b y -N t_CREB1_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 1 HEPG2 1 10 11 60 topRegions DNase

qsub -cwd -j y -q all.q -o ./m_CREB1_12.txt -b y -N m_CREB1_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 12 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_11.txt -b y -N m_CREB1_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 11 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_10.txt -b y -N m_CREB1_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 10 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_9.txt -b y -N m_CREB1_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 9 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_8.txt -b y -N m_CREB1_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 8 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_7.txt -b y -N m_CREB1_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 7 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_6.txt -b y -N m_CREB1_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 6 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_5.txt -b y -N m_CREB1_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 5 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_4.txt -b y -N m_CREB1_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 4 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_3.txt -b y -N m_CREB1_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 3 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_2.txt -b y -N m_CREB1_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 2 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_CREB1_1.txt -b y -N m_CREB1_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 1 HEPG2 1 10 500 550 middleRegions DNase

qsub -cwd -j y -q all.q -o ./b_CREB1_12.txt -b y -N b_CREB1_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 12 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_11.txt -b y -N b_CREB1_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 11 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_10.txt -b y -N b_CREB1_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 10 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_9.txt -b y -N b_CREB1_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 9 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_8.txt -b y -N b_CREB1_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 8 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_7.txt -b y -N b_CREB1_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 7 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_6.txt -b y -N b_CREB1_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 6 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_5.txt -b y -N b_CREB1_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 5 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_4.txt -b y -N b_CREB1_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 4 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_3.txt -b y -N b_CREB1_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 3 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_2.txt -b y -N b_CREB1_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 2 HEPG2 1 10 137000 137050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_CREB1_1.txt -b y -N b_CREB1_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R CREB1 1 HEPG2 1 10 137000 137050 bottomRegions DNase


########################## CHIPANALYSER DONE! ##########################

## PLOTTING ## create folders first.

# clean up folder (gzip *.fastq, gzip *.bdg, delete sam files)