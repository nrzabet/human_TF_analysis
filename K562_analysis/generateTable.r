rm(list=ls())

##get filenames
files <- read.csv('C:/Users/Romana/Dropbox/Documents/K562/SE/bdg_filenames_K562_ChIP_SE.txt',
                      sep ='\t', head = F, stringsAsFactors = F)
filenames <- data.frame(files[grep("treat_pileup.bdg", files$V1),])
rownames(filenames) <- NULL

##get narrowPeak filenames
narrow <- read.csv('C:/Users/Romana/Dropbox/Documents/K562/SE/K562NarrowPeakViable.txt',
                      sep ='\t', head = T, stringsAsFactors = F)
##get names of all TFs
TFs <- gsub("-human_K562_SE-peaks_treat_pileup.bdg", "", filenames[,1])

##col 1 - cell type will be pasted as the file name along with col 2
celltype <- rep("K562", nrow(filenames))

##col 2 - along with col 1, will form file name
outputName <- c()
for(i in TFs){
  outputName <- c(outputName, paste0(i, "_ChIP_PE"))
}

##col 3 - path to the PWM files
PFMs <- c()
for(j in TFs){
  PFMs <- c(PFMs, paste0("/home/rtpop/motifs/", j, paste0(".Rda")))
}

##col 4 - type of PWM - leave as matrix & convert all PWMs to the motifDb format
PFMtype <- rep("matrix", nrow(filenames))

##col 5 - path to bdg files
chip <- c()
for(k in filenames){
  chip <- c(chip, paste0("/home/rtpop/ENCODE_ChIP_data/K562/single_end/", k))
}

##col 6 - path to DNA sequence file
DNAsequence <- rep("/home/rtpop/hg38/hg38Sequence.Rda", nrow(filenames))

##col 7 - path to accessibility file
access <- rep("/home/rtpop/ENCODE_ChIP_data/K562/DNase-seq/dataset1/PE/quantile/K562_DNase_PE_broad_q01_quant0.1.Rda",
              nrow(filenames))

##col 8 - loci of interest - if NULL, peaks file should be provided
loci <- rep("NULL", nrow(filenames))

##col 9 -> no of regions to reduce to
reduce <- rep(60, nrow(filenames))

##col 10 path to narrowPeaks file -> only need if loci is NULL
peaks <- c()
for(l in narrow){
  peaks <- c(peaks, paste0("/home/rtpop/ENCODE_ChIP_data/K562/single_end/", l))
}

##col 11 - bin width for tiling the genome
lociWidth <- rep(20000, nrow(filenames))

##col 12 -> no of cores to use
cores <- rep(1, nrow(filenames))

##irrelevant, just keep as null; remove when you have time
col13 <- rep("NULL", nrow(filenames))

##col 14 - noise filtering method
noiseFil <- rep("sigmoid", nrow(filenames))

##col 15 - output directory - make sure it exists beforehand
outputDir <- rep("K562/single_end/ChIPanalyser/quant0.1/extraBM", nrow(filenames))

##col 16 - validation method - keep as AUC, the others are specified in
##ChIPanalPerformanalysys.r script
meth <- rep("AUC", nrow(filenames))

##as col 13
col17 <- rep("NULL", nrow(filenames))

## col 18 -> true if you want validation
validation <- rep("TRUE", nrow(filenames))

table <- cbind(celltype, outputName, PFMs, PFMtype, chip, DNAsequence, access,
              loci, reduce, peaks, lociWidth, cores, col13, noiseFil, outputDir,
              meth, col17, validation)

write.table(table, file = "K562_SE_quant0.1_extraBM.txt", sep = " ", col.names = FALSE,
            row.names = FALSE, quote = FALSE)
