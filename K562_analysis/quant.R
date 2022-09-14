## Loading libraries
library(rtracklayer)
library(GenomicRanges)

setwd('/home/rtpop/ENCODE_ChIP_data/K562/DNase-seq/dataset1/PE/')

### Loading DNA access .bdg
access <- import("K562_DNase_PE_broad_q01_treat_pileup.bdg", format="bedGraph")

# Set work directoty to wherever you want to save your files
setwd("/home/rtpop/ENCODE_ChIP_data/K562/DNase-seq/dataset1/PE/quantile")

# Creating Quantiled vector
quantVec<-c(seq(0,0.9, by=0.1), 0.95,0.99, 0.999)

# Looping over all quatiles
for(i in quantVec){
     # Your file name
     # Dont forget to change this for something more meaningful
     filename<-paste0("HepG2_DNase_PE_quant",i,".Rda")

     # subsetting Access based on quantiles
     # Make sure that the column with your enrichment score is called score
     # could be called scores.. can't remember now
     buffer<-reduce(access[which(access$score >= quantile(access$score,i))])
     save(buffer,file=filename)

}
