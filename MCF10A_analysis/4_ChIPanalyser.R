#!/usr/bin/env Rscript3.6.3
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Three arguments must be supplied", call.=FALSE)
} else {
  #transcription factor
  TF <- args[1]
  #Cell line
  cell_line <- args[2]
  #validation region
  regions <- args[3] #can be lost, new or maintained
}


setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/"))

library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPanalyser)


## Set parameters
noiseFil <- "sigmoid"
stepSize <- 100
chipSd <- 150
chipMean <- 150
boundMoleculesValues <- c(1, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 1e+06, 2e+06, 5e+06, 1e+07, 2e+07, 5e+07,1e+08)
lambdaValues <- seq(0.25,5, by=0.25)
PO <- parameterOptions(noiseFilter=noiseFil,chipSd=chipSd,chipMean=chipMean,stepSize=stepSize,boundMolecules=boundMoleculesValues,lambdaPWM=lambdaValues)



# Motifs from Romana's data
if(file.exists(paste0("/storage/st20d/HiC_Kc167/rtpop/motifs/",TF,".Rda"))){
  cat("loading /storage/st20d/HiC_Kc167/rtpop/motifs/",TF,".Rda file ...\n", sep ="")
  PFM <- get(load(paste0("/storage/st20d/HiC_Kc167/rtpop/motifs/",TF,".Rda")))
}else{
  stop(paste0("could not find motif for ",TF))
}


# hg38 objects
if(file.exists("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")){
  cat("loading /storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda file ...\n", sep ="")
  load("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")
} else{
  cat("recomputing DNAseqset object ...\n", sep ="")
  DNAseqset <- getSeq(BSgenome.Hsapiens.UCSC.hg38)[1:24]
  save(DNAseqset, file = "/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")
}


# accessibility data (we assume that everything is accessible and the TFs can bind in dense chromatin)
# access_data <- get(load(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/bdg/access_",cell_line,"_0.Rda")))


# ChIPseq data
regions_to_analyse <- get(load(paste0("peaks/",regions,"_regions_withNames.Rda")))


load(paste0("fpkm/FPKM_table_withOptimalParams.RData")) # FPKM table created in script 1
lambda <- fpkm_mat_withOptimalParams[TF,"optimal_lambda_K562"]
bm <- fpkm_mat_withOptimalParams[TF,"optimal_N_MCF10"]
PO <- parameterOptions(noiseFilter=noiseFil,chipSd=chipSd,chipMean=chipMean,stepSize=stepSize,boundMolecules=bm,lambdaPWM=lambda)

## Evaluation on the validation dataset
#control analysis
if(file.exists(paste0("rdaObjects/",TF,"/gw_ctrl_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/gw_ctrl_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/gw_ctrl_",TF,".Rda"))
} else {
  cat("recomputing gw_ctrl object ...\n", sep ="")
  GPP <- genomicProfiles(PFM=PFM,PFMFormat="matrix",BPFrequency=DNAseqset,PWMThreshold=0.7,lambdaPWM=lambda,boundMolecules=bm,stepSize=100,parameterOptions=PO)
  gw<-computeGenomeWideScores(GPP,DNAseqset,chromatinState=NULL,cores=1)
  save(gw, file=paste0("rdaObjects/",TF,"/gw_ctrl_",TF,".Rda"))
  rm(GPP)
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/pwm_ctrl_",regions,"Regions_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/pwm_ctrl_",regions,"Regions_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/pwm_ctrl_",regions,"Regions_",TF,".Rda"))
} else {
  cat("recomputing pwm_ctrl object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/gw_ctrl_",TF,".Rda")))
  pwm<-computePWMScore(genomicProfiles=gw,DNASequenceSet=DNAseqset,loci=regions_to_analyse,chromatinState=NULL,cores=1,parameterOptions=PO)
  rm(DNAseqset,gw,PO)
  save(pwm, file=paste0("rdaObjects/",TF,"/",regions,"/pwm_ctrl_",regions,"Regions_",TF,".Rda"))
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/occup_ctrl_",regions,"Regions_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/occup_ctrl_",regions,"Regions_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/occup_ctrl_",regions,"Regions_",TF,".Rda"))
} else {
  cat("recomputing occup_ctrl object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/",regions,"/pwm_ctrl_",regions,"Regions_",TF,".Rda")))
  occup<-computeOccupancy(pwm)
  rm(pwm)
  save(occup, file=paste0("rdaObjects/",TF,"/",regions,"/occup_ctrl_",regions,"Regions_",TF,".Rda"))
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_ctrl_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/predicted_chip_ctrl_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_ctrl_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))
} else {
  cat("recomputing predicted_chip_ctrl object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/",regions,"/occup_ctrl_",regions,"Regions_",TF,".Rda")))
  chip_ctrl<-computeChIPProfile(occup,regions_to_analyse,cores=10)
  rm(occup)
  save(chip_ctrl, file=paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_ctrl_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))
  rm(chip_ctrl)
}


## Set parameters
noiseFil <- "sigmoid"
stepSize <- 100
chipSd <- 150
chipMean <- 150
boundMoleculesValues <- c(1, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 1e+06, 2e+06, 5e+06, 1e+07, 2e+07, 5e+07,1e+08)
lambdaValues <- seq(0.25,5, by=0.25)
PO <- parameterOptions(noiseFilter=noiseFil,chipSd=chipSd,chipMean=chipMean,stepSize=stepSize,boundMolecules=boundMoleculesValues,lambdaPWM=lambdaValues)

load(paste0("fpkm/FPKM_table_withOptimalParams.RData"))
lambda <- fpkm_mat_withOptimalParams[TF,"optimal_lambda_K562"]
bm <- fpkm_mat_withOptimalParams[TF,"optimal_N_MCF10_HER2"]
PO <- parameterOptions(noiseFilter=noiseFil,chipSd=chipSd,chipMean=chipMean,stepSize=stepSize,boundMolecules=bm,lambdaPWM=lambda)

#Her2 overexpression analysis
if(file.exists(paste0("rdaObjects/",TF,"/gw_her2_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/gw_her2_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/gw_her2_",TF,".Rda"))
} else {
  cat("recomputing gw_her2 object ...\n", sep ="")
  load("/storage/st20d/HiC_Kc167/alessandra/hg38_objects/DNAseqset_hg38_gw.Rda")
  GPP <- genomicProfiles(PFM=PFM,PFMFormat="matrix",BPFrequency=DNAseqset,PWMThreshold=0.7,lambdaPWM=lambda,boundMolecules=bm,stepSize=100,parameterOptions=PO)
  gw<-computeGenomeWideScores(GPP,DNAseqset,chromatinState=NULL,cores=1)
  save(gw, file=paste0("rdaObjects/",TF,"/gw_her2_",TF,".Rda"))
  rm(GPP)
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/pwm_her2_",regions,"Regions_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/pwm_her2_",regions,"Regions_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/pwm_her2_",regions,"Regions_",TF,".Rda"))
} else {
  cat("recomputing pwm_her2 object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/gw_her2_",TF,".Rda")))
  pwm<-computePWMScore(genomicProfiles=gw,DNASequenceSet=DNAseqset,loci=regions_to_analyse,chromatinState=NULL,parameterOptions=PO,cores=1)
  rm(DNAseqset,gw,PO)
  save(pwm, file=paste0("rdaObjects/",TF,"/",regions,"/pwm_her2_",regions,"Regions_",TF,".Rda"))
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/occup_her2_",regions,"Regions_",TF,".Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/occup_her2_",regions,"Regions_",TF,".Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/occup_her2_",regions,"Regions_",TF,".Rda"))
} else {
  cat("recomputing occup_her2 object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/",regions,"/pwm_her2_",regions,"Regions_",TF,".Rda")))
  occup<-computeOccupancy(pwm)
  rm(pwm)
  save(occup, file=paste0("rdaObjects/",TF,"/",regions,"/occup_her2_",regions,"Regions_",TF,".Rda"))
}



if(file.exists(paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_her2_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))){
  cat("loading rdaObjects/",TF,"/",regions,"/predicted_chip_her2_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda file ...\n", sep ="")
  load(paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_her2_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))
} else {
  cat("recomputing predicted_chip_her2 object ...\n", sep ="")
  get(load(file=paste0("rdaObjects/",TF,"/",regions,"/occup_her2_",regions,"Regions_",TF,".Rda")))
  chip_her2<-computeChIPProfile(occup, regions_to_analyse, cores=1)
  rm(occup)
  save(chip_her2, file=paste0("rdaObjects/",TF,"/",regions,"/predicted_chip_her2_",regions,"Regions_",cell_line,"_",TF,"_gw.Rda"))
}
