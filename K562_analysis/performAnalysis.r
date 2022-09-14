##load relevant packages
library(ChIPanalyser)
library(GenomicRanges)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

######################################
########### Annex functions ##########
######################################

peaksLoading<-function(x){
    # Laoding files based on files extension

    if(grepl(x=x,pattern=".bed")& !grepl(x=x,pattern=".gff")){

          x<-read.table(x, stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }

    }else if(grepl(x=x, pattern=".Rda")){
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

performAnalysis <- function(TFList,ChIP,DNASequenceSet,
                          Access=Access,setSequence=setSequence,
                          reduce=reduce,peaks=peaks,tileSize=tileSize,cores=cores,
                          filename=filename,OP=OP,noiseFilter=noiseFil,method=methd){

#checking if the analysis was already run in case of multiple re-runs
if(file.exists(paste0(filename,"spearman_ChIPValidation.Rda")) == FALSE){
print("Loading DNA accessibility")
DNasePeaks <- peaksLoading(Access)
levels <- as.character(seqnames(DNasePeaks)) %in% c(paste0("chr",1:23),"chrX","chrY")
DNasePeaks <- GRanges(seqnames=as.character(seqnames(DNasePeaks))[levels],
ranges=IRanges(start(DNasePeaks[levels]), end(DNasePeaks[levels])))
rm(levels)

print("Loading genome sequence")
DNASequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38)

print("Loading PWM")
PFM <- get(load(TFList[[1]]))
GP <- genomicProfiles(PFM = PFM, PFMFormat= TFList[[2]],
                      BPFrequency = DNASequence)

number_of_regions_to_train <- 10

##set parameter options
PO <-parameterOptions(noiseFilter=noiseFilter,
                      chipSd=150,
                      chipMean=150,
                      lociWidth=tileSize,
                      stepSize=100)

if(file.exists(paste0(filename,"ChIPTraining.Rda"))){
  print("Processing ChIP")
  chipProfile <- get(load(paste0(filename,"ChIPTraining.Rda")))
}else{
print("Processing ChIP")
  chipProfile <- processingChIP(profile = ChIP, loci = setSequence, reduce = reduce,
                      peaks = peaks, chromatinState = DNasePeaks, cores = cores)
  save(chipProfile, file=paste0(filename,"ChIPTraining.Rda"))
}

##compute optimal
##get the training and validation regions    
if(length(loci(chipProfile))<=10){
  print("Fewer than 10 loci available")  
  TrainScore<-scores(chipProfile)
  ValidationScore<-scores(chipProfile)

  Trainloci<-loci(chipProfile)
  Validationloci<-loci(chipProfile)
}else if(length(loci(chipProfile))<60){
print("Fewer than 60 loci available")
print(paste0("Validation will be performed on", length(loci(chipProfile))))
TrainScore <- scores(chipProfile)[1:number_of_regions_to_train]
ValidationScore<-scores(chipProfile)[11:length(loci(chipProfile))]

Trainloci<-loci(chipProfile)[1:number_of_regions_to_train]
Validationloci <- loci(chipProfile)[11:length(loci(chipProfile))]
}else{
  TrainScore <- scores(chipProfile)[1:number_of_regions_to_train]
  ValidationScore<-scores(chipProfile)[11:60]

  Trainloci<-loci(chipProfile)[1:number_of_regions_to_train]
  Validationloci <- loci(chipProfile)[11:60]
}

#assigning the scores and loci for training regions
#a little bit hacky but faster than the proper way    
ChIPanalyser:::.scores(chipProfile)<-TrainScore
ChIPanalyser:::.loci(chipProfile)<-Trainloci

#getting rid of all the weird extra chromosomes    
chr <- as.character(seqnames(chipProfile@loci))
chr10 <- unique(chr)
DNASequence <- DNASequence[chr10]
rm(Trainloci)
rm(TrainScore)
save(Validationloci, file= paste0(filename, "ValidationLoci.Rda"))
save(ValidationScore, file= paste0(filename, "ValidationScore.Rda"))
rm(Validationloci)
rm(VatidationScore)

print("Computing optimal parameters")
optimal <- computeOptimal(genomicProfiles = GP,
                          DNASequenceSet = DNASequence,
                          ChIPScore = chipProfile,
                          chromatinState = DNasePeaks,
                          parameterOptions = PO,
                          cores = cores)

#get sequence for top 60 regions
if(!is.null(filename)){
save(optimal,file=paste0(filename,"OptimalOutputTraining.Rda"))
save(chipProfile, file=paste0(filename,"ChIPTraining.Rda"))
} else {
save(optimal,file=paste0("optimalOutputTraining.Rda"))
save(chipProfile, file=paste0("ChIPTraining.Rda"))
}

#get the validation regions
if(length(loci(chipProfile))<=10){
  ValidationScore<-scores(chipProfile)
  Validationloci<-loci(chipProfile)
}else{
ValidationScore<-scores(chipProfile)[11:60]
Validationloci <- loci(chipProfile)[11:60]
}

#validation with all the goodness of fit metrics
for(meth in method){

param<-optimal[[1]][[1]][[meth]]

lambda<-param[1]
bm<-param[2]
print("genomic profiles")
GPP<-genomicProfiles(PFM=PFM,PFMFormat=TFList[[2]],
                    BPFrequency=DNASequenceSet,PWMThreshold=0.7,lambdaPWM=lambda,
                    boundMolecules=bm,stepSize=100)
print("genome wide scores")
gw<-computeGenomeWideScores(GPP,DNASequence,DNasePeaks,cores=cores)
    
#assigning validation scores and loci
#again, a bit hacky but works    
ChIPanalyser:::.scores(chipProfile)<-ValidationScore
ChIPanalyser:::.loci(chipProfile)<-Validationloci

#getting rid of the weird chromosomes for these regions    
chr <- as.character(seqnames(chipProfile@loci))
chr10 <- unique(chr)
DNASequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38)[chr10]
    
print("pwm score")
pwm<-computePWMScore(gw,DNASequence,chipProfile,DNasePeaks,cores=cores)
occup<-computeOccupancy(pwm)

print("chip profile")
chip<-computeChIPProfile(occup,chipProfile,cores=cores)

print("accuracy estimate")
gof<-profileAccuracyEstimate(chip,chipProfile,cores=cores)
optimalList<-list("Occupancy"=occup,"ChIPProfile"=chip,"Gof"=gof)

if(!is.null(filename)){
    save(optimalList,file=paste0(filename,meth,"_","OptimalOutputValidation.Rda"))
    save(chipProfile, file=paste0(filename,meth,"_","ChIPValidation.Rda"))
} else {
    save(optimalList,file=paste0(meth,"_optimalOutputValidation.Rda"))
    save(chipProfile, file=paste0(meth,"_ChIPValidation.Rda"))
    }
  }
}else{
  print("File already exists")
}
}
