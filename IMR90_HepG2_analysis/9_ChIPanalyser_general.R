##############################################################################################
########################## General script for ChIPanalyser for hg38 ##########################
##############################################################################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=9) {
  stop("Nine arguments must be supplied", call.=FALSE)
} else {
  #transcription factor
  TF <- args[1]
  #DNA accessibility
  QDA_ID <- as.numeric(args[2])
  #Cell line
  cell_line <- args[3]
  #Regions to TRAIN, start and end regions
  regions_to_train_START <- as.numeric(args[4])
  regions_to_train_END <- as.numeric(args[5])
  #Regions to VALIDATE, start and end regions
  regions_to_test_START <- as.numeric(args[6])
  regions_to_test_END <- as.numeric(args[7])
  #specify range of validation regions 
  validation_regions <- args[8] #can either be topRegions, middleRegions or bottomRegions
  #DNA accessibility method
  DNAaccess <- args[9] #can be DNase, ATAC, MNase
}



setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/"))

library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
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



## Set parameters
tile_size <- 50000
lociWidth <- 50000
noiseFil <- "sigmoid"
stepSize <- 100
chipSd <- 150
chipMean <- 150
boundMoleculesValues <- c(1, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 
                          1e+05, 2e+05, 5e+05, 1e+06, 2e+06, 5e+06, 1e+07, 2e+07, 5e+07,1e+08)
lambdaValues <- seq(0.25,5, by=0.25)
PO <- parameterOptions(noiseFilter=noiseFil,
                      chipSd=chipSd,
                      chipMean=chipMean,
                      lociWidth=lociWidth,
                      stepSize=stepSize,
                      boundMolecules=boundMoleculesValues,
                      lambdaPWM = lambdaValues)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95, 0.99)
regions_to_select <- max(c(regions_to_train_START, regions_to_train_END,
                           regions_to_test_START, regions_to_test_END)) # the max of the 4 values



## hg38 objects
if(file.exists("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda file ...\n", sep ="")
  load("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")
} else{
  cat("recomputing DNAseqset object ...\n", sep ="")
  DNAseqset <- getSeq(BSgenome.Hsapiens.UCSC.hg38)[1:24]
  save(DNAseqset, file = "/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")
}

if(file.exists("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangegoodbeforetiles_hg38_gw.Rda")){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangegoodbeforetiles_hg38_gw.Rda file ...\n", sep ="")
  load("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangegoodbeforetiles_hg38_gw.Rda")
} else{
  cat("recomputing grangegood object ...\n", sep ="")
  grangegood <- GRanges(seqnames=names(BSgenome.Hsapiens.UCSC.hg38)[1:24],
                        ranges=IRanges(start=1, end=width(getSeq(BSgenome.Hsapiens.UCSC.hg38))[1:24]))
  save(grangegood, file = "/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangegoodbeforetiles_hg38_gw.Rda")
}

if(file.exists("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangesobjectaftertiles_hg38_gw.Rda")){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangesobjectaftertiles_hg38_gw.Rda file ...\n", sep ="")
  load("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangesobjectaftertiles_hg38_gw.Rda")
} else{
  cat("recomputing grangesobjecttiles object ...\n", sep ="")
  grangesobjecttiles <- unlist(tile(GRanges(seqnames=names(BSgenome.Hsapiens.UCSC.hg38)[1:24],
                                            ranges=IRanges(start=1, end=width(getSeq(BSgenome.Hsapiens.UCSC.hg38))[1:24])),
                                    width=tile_size))
  save(grangesobjecttiles, file="/storage/st20d/HiC_Kc167/alessandra/hg38_objects/grangesobjectaftertiles_hg38_gw.Rda")
}



## ChIPseq data
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))
} else{
  cat("recomputing profile_data object ...\n", sep ="")
  profile_data <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/",TF,"-human_",cell_line,"_SE_trim-peaks_treat_pileup.bdg"), format = "bedGraph")
  save(profile_data,file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))
}

if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_peaks_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_peaks_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_peaks_gw.Rda"))
} else{
  cat("recomputing peaks_data object ...\n", sep ="")
  peaks_data <- peaksLoading(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/",TF,"-human_",cell_line,"_SE_trim-peaks_peaks.narrowPeak"))
  save(peaks_data, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_peaks_gw.Rda"))
}



## Motifs (already have those from 4-prechipanalyser.R)
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",TF,".Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",TF,".Rda file ...\n", sep ="")
  PFM <- load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",TF,".Rda"))
} else{
  stop(paste0("could not find motif for ",TF))
}

load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/motifs/",TF,".Rda"))



## Accessibility data for specific QDA_ID
access_data<-paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/",DNAaccess,"/access_",cell_line,"_",quantVec[QDA_ID],".Rda")
if(file.exists(access_data)){
  cat("recomputing access_data object ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/",DNAaccess,"/access_",cell_line,"_",quantVec[QDA_ID],".Rda"))
} else{
  stop(paste0("cannot load accessibility data: ",access_data))
}



## ChIPanalyser
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_MM_chipscore_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing chipscore object ...\n", sep ="")
  chipscore <-processingChIP(profile=profile_data,
                             loci=grangesobjecttiles,
                             reduce=regions_to_select,
                             peaks= peaks_data,
                             cores=1,
                             chromatinState=access_data,
                             parameterOptions=PO)
  save(chipscore, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))
}
rm(profile_data, grangesobjecttiles)

if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/analysisRegions_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/analysisRegions_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/analysisRegions_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  analysisRegions <- chipscore 
  TrainScore <- scores(chipscore)[regions_to_train_START:regions_to_train_END]
  ValidationScore<-scores(chipscore)[regions_to_test_START:regions_to_test_END]
  TrainLoci<-loci(chipscore)[regions_to_train_START:regions_to_train_END] 
  ValidationLoci <- loci(chipscore)[regions_to_test_START:regions_to_test_END]
  ChIPanalyser:::.scores(analysisRegions)<-TrainScore
  ChIPanalyser:::.loci(analysisRegions)<-TrainLoci # Compute optimal
  save(analysisRegions, TrainScore, TrainLoci, ValidationScore, ValidationLoci, 
       file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/analysisRegions_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}
rm(chipscore)



## PWM scores
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_GP_access_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_GP_access_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_GP_access_gw.Rda"))
} else{
  cat("recomputing TF_GP object ...\n", sep ="")
  TF_GP<-genomicProfiles(PFM=PFM,
                         PFMFormat="raw",
                         BPFrequency=DNAseqset,
                         ChIPScore=NULL,
                         boundMolecules=boundMoleculesValues,
                         lambdaPWM=lambdaValues)
  save(TF_GP, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_GP_access_gw.Rda"))
}



## Keep chromosomes to analyse and remove computationally-heavy objects
y<- unique(analysisRegions@loci@seqnames)				
DNAseqset_sub<-getSeq(BSgenome.Hsapiens.UCSC.hg38)[y]
rm(DNAseqset)
rm(peaksLoading)
rm(TrainScore, TrainLoci, ValidationScore, ValidationLoci)
rm(peaks_data)



## Optimal parameters
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing optimal_TF object ...\n", sep ="")
  optimal_TF<-suppressWarnings(computeOptimal(genomicProfiles=TF_GP,
                                              DNASequenceSet=DNAseqset_sub, 
                                              ChIPScore=analysisRegions,
                                              chromatinState=access_data,
                                              cores=1))
  save (optimal_TF, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}



## Evaluation on the validation dataset
if(file.exists(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing optimal_validation_TF object ...\n", sep ="")
  param <- optimal_TF[[1]][[1]][["MSE"]]
  lambda <- param[1]
  bm <- param[2] # compute global scores
  GPP <- genomicProfiles(PFM=PFM,PFMFormat="matrix",
                         BPFrequency=DNAseqset_sub,PWMThreshold=0.7,
                         lambdaPWM=lambda,boundMolecules=bm,stepSize=100,
                         parameterOptions=PO)
  rm(optimal_TF, analysisRegions)
  
  if(quantVec[QDA_ID] > 0){
    gw<-computeGenomeWideScores(GPP,DNAseqset_sub,chromatinState=access_data,cores=1) 
  } else{
    gw<-computeGenomeWideScores(GPP,DNAseqset_sub,chromatinState=NULL,cores=1)
  }
  
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/analysisRegions_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda")) #load TF-specific analysisRegion object
  
  ChIPanalyser:::.scores(analysisRegions)<-ValidationScore
  ChIPanalyser:::.loci(analysisRegions)<-ValidationLoci
  
  cat("recomputing optimal_validation_TF object ...\n", sep ="")
  
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_peaks_gw.Rda")) #load TF-specific peaks object
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")) # load general (non TF-specific) hg38 object
  rm(DNAseqset_sub)
  load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda")) #load TF-specific peaks object
  pwm<-computePWMScore(genomicProfiles=gw,DNASequenceSet=DNAseqset, 
                       loci=analysisRegions,chromatinState=access_data,
                       parameterOptions=PO,cores=1)
  rm(DNAseqset, access_data, gw)
  occup<-computeOccupancy(pwm)
  chip<-computeChIPProfile(occup, analysisRegions, cores=1)
  gof<-profileAccuracyEstimate(chip, analysisRegions, cores=1)
  optimal_validation_TF<-list("Occupancy"=occup,"ChIPProfile"=chip,"GOF"=gof)
  save(optimal_validation_TF, file=paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ChIPanalyser/rdaObjects/",DNAaccess,"/",TF,"/",validation_regions,"/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}