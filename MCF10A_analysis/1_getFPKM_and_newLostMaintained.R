############# Get FPKMs for TFs of interest in MCF10 and K562 ############# 

# To get the files go to:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1420579 (MCF10)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1940175 (K562)
# Hit download on "http" at the bottom of the page. Move these files to alessandra/MCF10 folder
# Unzip them with gunzip filename.txt.gz
# R3.6.0

setwd("/storage/st20d/HiC_Kc167/alessandra/MCF10/fpkm/")

TFs <- c("ATF1","ETV6","MYC","SRF","JUN","JUND","GTF3C2","NFATC3","CHD7","SOX6","TRIM28","YBX1")
# TFs I couldn't find: ZBTB,CHAMP1,ZNF,ZBT
cellLines <- c("MCF10_fpkm","K562_fpkm")

MCF10 <- read.table("MCF10_FPKM.txt", stringsAsFactors=F, header=T)
K562 <- read.table("K562_FPKM.txt", stringsAsFactors=F, header=T)

fpkm_mat <- matrix(data=NA, nrow=12, ncol=2)
rownames(fpkm_mat) <- TFs
colnames(fpkm_mat) <- cellLines

for(i in TFs){
  fpkm <- MCF10[which(MCF10$tracking_id == i),"FPKM"]
  fpkm_mat[i,1] <- fpkm
}
for(i in TFs){
  fpkm <- K562[which(K562$tracking_id == i),"FPKM"]
  fpkm_mat[i,2] <- fpkm
}

save(fpkm_mat, file="FPKM_table.RData")


# Get optimal params for K562 from Romana's data for these TFs
library(ChIPanalyser)
TFs <- c("ATF1","ETV6","MYC","SRF","JUN","JUND","GTF3C2","NFATC3","CHD7","SOX6","TRIM28","YBX1")

romanaData <- matrix(data=NA, nrow=12, ncol=2)
rownames(romanaData) <- TFs
colnames(romanaData) <- c("optimal_lambda_K562","optimal_N_K562")

for(i in TFs){
  if(file.exists(paste0("/storage/st20d/HiC_Kc167/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant0/",i,"_ChIP_SE_K562_top10_OptimalOutputTraining.Rda"))){ # paired ends
    setwd("/storage/st20d/HiC_Kc167/rtpop/ENCODE_ChIP_data/K562/single_end/ChIPanalyser/quant0/")
    load(paste0(i,"_ChIP_SE_K562_top10_OptimalOutputTraining.Rda"))
    romanaData[i,] <- optimal$Optimal$OptimalParameters$MSE
  }else{
    setwd("/storage/st20d/HiC_Kc167/rtpop/ENCODE_ChIP_data/K562/pair_end/ChIPanalyser/quant0/") # single ends
    load(paste0(i,"_ChIP_PE_K562_top10_OptimalOutputTraining.Rda"))
    romanaData[i,] <- optimal$Optimal$OptimalParameters$MSE
  }}

load("/storage/st20d/HiC_Kc167/alessandra/MCF10/fpkm/FPKM_table.RData") 
fpkm_mat_withOptimalParams <- cbind(fpkm_mat,romanaData)


# compute optimal N for MCF10
optimal_N_MCF10 <- (fpkm_mat_withOptimalParams[,4] * (fpkm_mat_withOptimalParams[,1]/fpkm_mat_withOptimalParams[,2]))
fpkm_mat_withOptimalParams <- cbind(fpkm_mat_withOptimalParams,optimal_N_MCF10)


# compute optimal N for MCF10 and HER2 by rescaling the optimal_N_MCF10 based on the numbers Ateeq sent me
fpkm_mat_withOptimalParams <- rbind(fpkm_mat_withOptimalParams,fpkm_mat_withOptimalParams[8,]) # because NFATC3 has 2 values in Ateeqs table
fpkm_mat_withOptimalParams <- rbind(fpkm_mat_withOptimalParams,fpkm_mat_withOptimalParams[11,]) # because TRIM28 has 3 values in Ateeqs table
fpkm_mat_withOptimalParams <- rbind(fpkm_mat_withOptimalParams,fpkm_mat_withOptimalParams[11,])

rownames(fpkm_mat_withOptimalParams) <- c("ATF1","ETV6","MYC","SRF","JUN","JUND","GTF3C2","NFATC3","CHD7","SOX6","TRIM28","YBX1","NFATC3_2","TRIM28_2","TRIM28_3")
fpkm_mat_withOptimalParams <- fpkm_mat_withOptimalParams[order(rownames(fpkm_mat_withOptimalParams)),] # sort matri so that TFs are in alphabetical order

sevenHours <- c(-0.69269532723869, -0.497177017170269, -1.01011800214394, 0.662965012722427, 0.184942037083036,
                -0.520545013192049, -0.946145660801075, -0.962257906008448, 1.00734984986277, -0.4589208756448115,
                1.32701058635314, -1.67434942900012, 0.524829018480776, -0.553073064136331, 0.570315724756754) # these come from Ateeqs excel table column K. Make sure they are in order based on the matrix's rownames (the TFs).

optimal_N_MCF10_HER2 <- fpkm_mat_withOptimalParams[,5] * (2 ^ sevenHours) # optimal_N_MCF10 x (2 to the power of column K i.e. sevenHours)

fpkm_mat_withOptimalParams <- cbind(fpkm_mat_withOptimalParams,optimal_N_MCF10_HER2)
save(fpkm_mat_withOptimalParams, file="/storage/st20d/HiC_Kc167/alessandra/MCF10/fpkm/FPKM_table_withOptimalParams.RData")