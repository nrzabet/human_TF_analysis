#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#############################################################################################
### What to analyse
#############################################################################################


if (length(args)!=3) {
  stop("Three arguments must be supplied", call.=FALSE)
} else {

  #transcription factor
  TF <- args[1]

  #DNA accessibility
  QDA_ID <- as.numeric(args[2])

  #Cell line
  cell_line <- args[3]
}

################################################################################
# libraries and functions
################################################################################
setwd("/storage/projects/ZabetLab/mm10_ENCODE/scripts/")

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

################################################################################
# parameters
################################################################################
tile_size <- 50000
regions_to_select <- 60
number_of_regions_to_train <- 10
lociWidth <- 50000
noiseFil <- "sigmoid"
stepSize <- 100
chipSd <- 150
chipMean <- 150
boundMoleculesValues <- c(1, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 
						1e+05, 2e+05, 5e+05, 1e+06, 2e+06, 5e+06, 1e+07, 2e+07, 5e+07,1e+08)
lambdaValues <- seq(0.25,5, by=0.25)
PO <-parameterOptions(noiseFilter=noiseFil,
                      chipSd=chipSd,
                      chipMean=chipMean,
                      lociWidth=lociWidth,
                      stepSize=stepSize,
                      boundMolecules=boundMoleculesValues,
                      lambdaPWM = lambdaValues)

quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)

################################################################################
# mm10 objects
################################################################################
if(file.exists("../objects/DNAseqset_mm10_gw.Rda")){
  cat("loading ../objects/DNAseqset_mm10_gw.Rda file ...\n", sep ="")
  load("../objects/DNAseqset_mm10_gw.Rda")
} else{
  cat("recomputing DNAseqset object ...\n", sep ="")
  DNAseqset <- getSeq(BSgenome.Mmusculus.UCSC.mm10)[1:21]
	save(DNAseqset, file = "../objects/DNAseqset_mm10_gw.Rda")
}

if(file.exists("../objects/grangegoodbeforetiles_mm10_gw.Rda")){
  cat("loading ../objects/grangegoodbeforetiles_mm10_gw.Rda file ...\n", sep ="")
  load("../objects/grangegoodbeforetiles_mm10_gw.Rda")
} else{
  cat("recomputing grangegood object ...\n", sep ="")
	grangegood <- GRanges(seqnames=names(BSgenome.Mmusculus.UCSC.mm10)[1:21],
                        ranges=IRanges(start=1, end=width(getSeq(BSgenome.Mmusculus.UCSC.mm10))[1:21]))
	save(grangegood, file = "../objects/grangegoodbeforetiles_mm10_gw.Rda")
}

if(file.exists("../objects/grangesobjectaftertiles_mm10_gw.Rda")){
  cat("loading ../objects/grangesobjectaftertiles_mm10_gw.Rda file ...\n", sep ="")
  load("../objects/grangesobjectaftertiles_mm10_gw.Rda")
} else{
	cat("recomputing grangesobjecttiles object ...\n", sep ="")
  grangesobjecttiles <- unlist(tile(GRanges(seqnames=names(BSgenome.Mmusculus.UCSC.mm10)[1:21],
                      ranges=IRanges(start=1, end=width(getSeq(BSgenome.Mmusculus.UCSC.mm10))[1:21])),
              width=tile_size))
	save(grangesobjecttiles, file="../objects/grangesobjectaftertiles_mm10_gw.Rda")
}

################################################################################
# ChIPseq data
################################################################################
if(file.exists(paste0("../objects/set2/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))){
  cat("loading ../objects/set2/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))
} else{
  cat("recomputing profile_data object ...\n", sep ="")
	profile_data <- import(paste0("../data/set2/ChIPseq_",cell_line,"_",TF,"_peaks_treat_pileup.bdg"), format = "bedGraph")
  save(profile_data,file=paste0("../objects/set2/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))
}


if(file.exists(paste0("../objects/set2/",cell_line,"_",TF,"_peaks_gw.Rda"))){
  cat("loading ../objects/set2/",cell_line,"_",TF,"_peaks_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/",cell_line,"_",TF,"_peaks_gw.Rda"))
} else{
  cat("recomputing peaks_data object ...\n", sep ="")
	peaks_data <- peaksLoading(paste0("../data/set2/ChIPseq_",cell_line,"_",TF,"_peaks_peaks.narrowPeak"))
	save(peaks_data, file=paste0("../objects/set2/",cell_line,"_",TF,"_peaks_gw.Rda"))
}

if(file.exists(paste0("../objects/set2/motifs/",TF,".Rda"))){
  cat("loading ../objects/set2/motifs/",TF,".Rda file ...\n", sep ="")
  load(paste0("../objects/set2/motifs/",TF,".Rda"))
} else{
  cat("recomputing PFM object ...\n", sep ="")
	PFM<-query(query(query(MotifDb,"mmusculus"),paste0(TF)) ,"HOCOMOCOv10")
               PFM <- PFM[[1]]
	save(PFM, file = paste0("../objects/set2/motifs/",TF,".Rda"))
}

################################################################################
# accessibility data
################################################################################
access_filename<-paste0("../objects/set2/access_",cell_line,"_",quantVec[QDA_ID],"_gw.Rda")



if(file.exists(access_filename)){
  cat("recomputing access_data object ...\n", sep ="")
	access_data <- peaksLoading(access_filename)
} else{
  stop(paste0("cannot load accessibility data: ",access_filename))
}



################################################################################
# ChIPanalyser
################################################################################
# add cell line
if(file.exists(paste0("../objects/set2/",cell_line,"_",TF,"MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading ../objects/set2/",cell_line,"_",TF,"MM_chipscore_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/",cell_line,"_",TF,"MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing chipscore object ...\n", sep ="")
  chipscore <-processingChIP(profile=profile_data,
               loci=grangesobjecttiles,
               reduce=regions_to_select,
               peaks= peaks_data ,
               cores=1,
			   chromatinState=access_data,
			   parameterOptions=PO)
  save(chipscore, file = paste0("../objects/set2/",cell_line,"_",TF,"MM_chipscore_",quantVec[QDA_ID],"_gw.Rda"))
}
rm(profile_data, grangesobjecttiles)


if(file.exists(paste0("../objects/set2/top50training",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading ../objects/set2/top50training",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/top50training",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  ##compute optimal
  top50training <- chipscore

  TrainScore <- scores(chipscore)[1:number_of_regions_to_train]
  ValidationScore<-scores(chipscore)[(number_of_regions_to_train+1):regions_to_select]

  TrainLoci<-loci(chipscore)[1:number_of_regions_to_train]
  ValidationLoci <- loci(chipscore)[(number_of_regions_to_train+1):regions_to_select]

  ChIPanalyser:::.scores(top50training)<-TrainScore
  ChIPanalyser:::.loci(top50training)<-TrainLoci
  save(top50training, TrainScore, TrainLoci, ValidationScore, ValidationLoci, 
	file=paste0("../objects/set2/top50training",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}
rm(chipscore)



# PWM scores
if(file.exists(paste0("../objects/set2/",cell_line,"_",TF,"_GP_access_gw.Rda"))){
  cat("loading ../objects/set2/",cell_line,"_",TF,"_GP_access_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/",cell_line,"_",TF,"_GP_access_gw.Rda"))
} else{
  cat("recomputing TF_GP object ...\n", sep ="")
	TF_GP<-genomicProfiles(PFM=PFM,
                           PFMFormat="raw",
                           BPFrequency=DNAseqset,
                           ChIPScore=top50training,
						   boundMolecules=boundMoleculesValues,
	                       lambdaPWM=lambdaValues)
	save(TF_GP, file=paste0("../objects/set2/",cell_line,"_",TF,"_GP_access_gw.Rda"))
}

y<- unique(top50training@loci@seqnames)				
						
DNAseqset_sub<-getSeq(BSgenome.Mmusculus.UCSC.mm10)[y]
  rm(DNAseqset)
 
						
rm(peaksLoading)
  rm(TrainScore, TrainLoci, ValidationScore, ValidationLoci)
  rm(peaks_data)

#optimal paramters
if(file.exists(paste0("../objects/set2/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading ../objects/set2/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing optimal_TF object ...\n", sep ="")
  optimal_TF<-suppressWarnings(computeOptimal(genomicProfiles=TF_GP,
                                                     DNASequenceSet=DNAseqset_sub, 
                                                     ChIPScore=top50training,
                                                     chromatinState=access_data,
													 cores=1))
  save (optimal_TF, file=paste0("../objects/set2/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}


#evaluation on the validation dataset
if(file.exists(paste0("../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
  cat("loading ../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda file ...\n", sep ="")
  load(paste0("../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
} else{
  cat("recomputing optimal_validation_TF object ...\n", sep ="")

  ## compute global scores
  param <- optimal_TF[[1]][[1]][["MSE"]]
  lambda <- param[1]
  bm <- param[2]


  GPP <- genomicProfiles(PFM = PFM , PFMFormat = "matrix",
                        BPFrequency = DNAseqset_sub, PWMThreshold=0.7,
                        lambdaPWM=lambda,boundMolecules=bm,stepSize=100,
                        parameterOptions=PO)

  objects <- sort( sapply(ls(),function(x){object.size(get(x))}))
  print(objects)

rm(optimal_TF, top50training)

  if(quantVec[QDA_ID] > 0){
    gw<-computeGenomeWideScores(GPP, DNAseqset_sub, chromatinState = access_data, cores=1) 
  } else{
    gw<-computeGenomeWideScores(GPP, DNAseqset_sub, chromatinState = NULL, cores=1)
  }
  
  load(paste0("../objects/set2/top50training",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
  

  ChIPanalyser:::.scores(top50training)<-ValidationScore
  ChIPanalyser:::.loci(top50training)<-ValidationLoci

  cat("recomputing optimal_validation_TF object ...\n", sep ="")

  load(paste0("../objects/set2/",cell_line,"_",TF,"_peaks_gw.Rda"))
   load("../objects/DNAseqset_mm10_gw.Rda")
rm(DNAseqset_sub)
  load(paste0("../objects/set2/",cell_line,"_",TF,"_ChIPseqprofile_gw.Rda"))
  pwm<-computePWMScore(genomicProfiles = gw, DNASequenceSet = DNAseqset, 
                       loci = top50training, chromatinState = access_data,
                       parameterOptions = PO, cores=1)
  
 rm(DNAseqset, access_data, gw)
  occup<-computeOccupancy(pwm)
  

  chip<-computeChIPProfile(occup, top50training, cores=1)

  gof<-profileAccuracyEstimate(chip, top50training, cores=1)

  optimal_validation_TF<-list("Occupancy"=occup,"ChIPProfile"=chip,"GOF"=gof)

  save(optimal_validation_TF, file=paste0("../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
}


optimal_validation_TF$GOF@profiles[[1]][[1]]["AUCMean"]
