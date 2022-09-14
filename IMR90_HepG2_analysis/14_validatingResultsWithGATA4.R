######### Validating results with known pioneer TFs in HepG2 with ChIPanalyser #########
## GATA4
## https://www.encodeproject.org/matrix/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=TF+ChIP-seq&target.label=GATA4

## PREPROCESSING ##
## Open link above and select fastq only. Hit download. Open "files.txt", run first line on Google. It will download a file called "metadata".
## Change metadata name to match cell line e.g. metadata_cell_line_TF.tsv
## Move everything to cluster. Open screen, qrsh, go to right directory
xargs -L 1 curl -O -J -L < files.txt 
gunzip *.gz # if the files are zipped

## In the same working directory type:
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
## Now go inside Trimmomatic-0.39 and move the whole content outside of its folder, into preProcessing folder.

R3.5.0
cell_line <- "HEPG2"
srr1 <- "ENCFF193ROJ"
srr2 <- "ENCFF163SRP" # fastq names from files

commands1<- c("#!/bin/bash\n\n")
commands1 <- paste0("cat ",srr1,".fastq ", srr2, ".fastq > GATA4_",cell_line,".fastq", collapse =  "\n")
write(commands1, file="download_GATA4_files_1.sh")

commands2 <- c("#!/bin/bash\n\n")
fastqc<-paste0("fastqc GATA4_",cell_line,".fastq")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U GATA4_",cell_line,".fastq -S 'GATA4_",cell_line,"_chip.sam' > 'GATA4_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t GATA4_",cell_line,"_chip.sam -f SAM -g 2.9e9 -n GATA4_",cell_line,"_peaks -B -q 0.05 > 'GATA4_",cell_line, "_peaksLog.log'")
commands2<-paste0(fastqc, "\n", bowtie, "\n", macs2, "\n")
write(commands2, file="download_GATA4_files_2.sh")

q()
n

chmod 770 download_GATA4_files_1.sh
chmod 770 download_GATA4_files_2.sh
./download_GATA4_files_1.sh
./download_GATA4_files_2.sh

########################## PRE-PROCESSING DONE! ##########################  



## FINDING NUMBER OF CHIPSEQ PEAKS ##
R3.5.0
cell_line <- "HepG2"
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

tfs <- "GATA4"
cell_line <- "HEPG2"

peaksPlot<- function(tfs){
  peaks_data <- peaksLoading(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/",
                                    tfs,"_",cell_line,"_peaks_peaks.narrowPeak"))
  x<- length(peaks_data)
}
peaks <- sapply(tfs, peaksPlot) # number of peaks for GATA4: 36300


########################## NUMBER OF CHIP-SEQ DONE! ##########################



## GET MOTIFS ##
## screen, qrsh, cd /, go to IMR90/preProcessing folder
ls GATA4_HEPG2_peaks_treat_pileup.bdg > HEPG2_GATA4_bdg.txt
mv /storage/st20d/HiC_Kc167/alessandra/HEPG2/preProcessing/HEPG2_GATA4_bdg.txt /storage/st20d/HiC_Kc167/alessandra/HEPG2/ChIPanalyser/motifs/
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/ChIPanalyser/motifs/
R3.5.0

library(MotifDb)
cell_line <- "HEPG2"
TFnames <- "GATA4"

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
qsub -cwd -j y -q all.q -o ./t_GATA4_12.txt -b y -N t_GATA4_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 12 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_11.txt -b y -N t_GATA4_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 11 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_10.txt -b y -N t_GATA4_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 10 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_9.txt -b y -N t_GATA4_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 9 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_8.txt -b y -N t_GATA4_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 8 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_7.txt -b y -N t_GATA4_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 7 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_6.txt -b y -N t_GATA4_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 6 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_5.txt -b y -N t_GATA4_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 5 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_4.txt -b y -N t_GATA4_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 4 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_3.txt -b y -N t_GATA4_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 3 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_2.txt -b y -N t_GATA4_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 2 HEPG2 1 10 11 60 topRegions DNase
qsub -cwd -j y -q all.q -o ./t_GATA4_1.txt -b y -N t_GATA4_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 1 HEPG2 1 10 11 60 topRegions DNase

qsub -cwd -j y -q all.q -o ./m_GATA4_12.txt -b y -N m_GATA4_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 12 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_11.txt -b y -N m_GATA4_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 11 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_10.txt -b y -N m_GATA4_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 10 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_9.txt -b y -N m_GATA4_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 9 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_8.txt -b y -N m_GATA4_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 8 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_7.txt -b y -N m_GATA4_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 7 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_6.txt -b y -N m_GATA4_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 6 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_5.txt -b y -N m_GATA4_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 5 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_4.txt -b y -N m_GATA4_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 4 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_3.txt -b y -N m_GATA4_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 3 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_2.txt -b y -N m_GATA4_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 2 HEPG2 1 10 500 550 middleRegions DNase
qsub -cwd -j y -q all.q -o ./m_GATA4_1.txt -b y -N m_GATA4_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 1 HEPG2 1 10 500 550 middleRegions DNase

qsub -cwd -j y -q all.q -o ./b_GATA4_12.txt -b y -N b_GATA4_12 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 12 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_11.txt -b y -N b_GATA4_11 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 11 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_10.txt -b y -N b_GATA4_10 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 10 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_9.txt -b y -N b_GATA4_9 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 9 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_8.txt -b y -N b_GATA4_8 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 8 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_7.txt -b y -N b_GATA4_7 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 7 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_6.txt -b y -N b_GATA4_6 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 6 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_5.txt -b y -N b_GATA4_5 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 5 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_4.txt -b y -N b_GATA4_4 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 4 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_3.txt -b y -N b_GATA4_3 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 3 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_2.txt -b y -N b_GATA4_2 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 2 HEPG2 1 10 22000 22050 bottomRegions DNase
qsub -cwd -j y -q all.q -o ./b_GATA4_1.txt -b y -N b_GATA4_1 Rscript ./ChIPanalyser_general_forKnownPioneers.R GATA4 1 HEPG2 1 10 22000 22050 bottomRegions DNase

########################## CHIPANALYSER DONE! ##########################

## PLOTTING ## create folders first.

# clean up folder (gzip *.fastq, gzip *.bdg, delete sam files)