##################### HepG2 DNA access preprocessing #####################


######## ATAC #######
# ATAC-seq https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5214186
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/ATAC
fastq-dump SRR14103581
R3.5.0
cell_line <- "HEPG2" 
srr1 <- "SRR14103581"

commands1 <- c("#!/bin/bash\n\n")
cat <- paste0("cat ",srr1,".fastq > ATAC_",cell_line,".fastq", collapse =  "\n")
fastqc<-paste0("fastqc ATAC_",cell_line,".fastq -o fastqc_raw/")
cutadapt<-paste0("cutadapt -u 19 -o ATAC_",cell_line,"_forward.fq.gz" , " ATAC_",cell_line, ".fastq")
fastqc2<-paste0("fastqc ATAC_",cell_line,"_forward.fq.gz -o fastqc_raw/") # continue here, run this.
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U ATAC_",cell_line,"_forward.fq.gz -S 'ATAC_",cell_line,"_chip.sam' > 'ATAC_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t ATAC_",cell_line,"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n ATAC_",cell_line,"_peaks -B --nomodel > 'ATAC_",cell_line, "_peaksLog.log'")
commands1<-paste0(cat, "\n", fastqc, "\n", cutadapt, "\n", fastqc2, "\n", bowtie, "\n", macs2, "\n")
write(commands1, file="download_ATAC_files_1.sh")

quit()
n
chmod 770 download_ATAC_files_1.sh
./download_ATAC_files_1.sh

# cd /, R3.6.0
cell_line <- "HEPG2"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC/ATAC_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  save(access_data,file=access_cell_line)
}





######## DNASE #######
# DNase-seq https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2400286
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/DNase
fastq-dump SRR5048293
fastq-dump SRR5048294
fastq-dump SRR5048295
R3.5.0
cell_line <- "HEPG2"
srr1 <- "SRR5048293"
srr2 <- "SRR5048294"
srr3 <- "SRR5048295"

commands1 <- c("#!/bin/bash\n\n")
cat <- paste0("cat ",srr1,".fastq ",srr2,".fastq ",srr3,".fastq"," > DNase_",cell_line,".fastq", collapse =  "\n")
fastqc<-paste0("fastqc DNase_",cell_line,".fastq -o fastqc_raw/")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U DNase_",cell_line,".fastq -S 'DNase_",cell_line,"_chip.sam' > 'DNase_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t DNase_",cell_line,"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n DNase_",cell_line,"_peaks -B --nomodel > 'DNase_",cell_line, "_peaksLog.log'")
commands1<-paste0(cat, "\n", fastqc, "\n", bowtie, "\n", macs2, "\n")
write(commands1, file="download_DNase_files_1.sh")

quit()
n
chmod 770 download_DNase_files_1.sh
./download_DNase_files_1.sh

# cd /, R3.6.0
cell_line <- "HEPG2"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase/DNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  save(access_data,file=access_cell_line)
}






####### MNASE ######
# MNase-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3718063
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/MNase
fastq-dump SRR5048293
fastq-dump SRR5048294
fastq-dump SRR5048295
R3.5.0
cell_line <- "HEPG2"
srr1 <- "SRR8882563"
srr2 <- "SRR8882564"
srr3 <- "SRR8882565"
srr4 <- "SRR8882566"

commands1 <- c("#!/bin/bash\n\n")
cat <- paste0("cat ",srr1,".fastq ",srr2,".fastq ",srr3,".fastq ",srr4,".fastq"," > MNase_",cell_line,".fastq", collapse =  "\n")
fastqc<-paste0("fastqc MNase_",cell_line,".fastq -o fastqc_raw/")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U MNase_",cell_line,".fastq -S 'MNase_",cell_line,"_chip.sam' > 'MNase_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t MNase_",cell_line,"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n MNase_",cell_line,"_peaks -B --nomodel > 'MNase_",cell_line, "_peaksLog.log'")
commands1<-paste0(cat, "\n", fastqc, "\n", bowtie, "\n", macs2, "\n")
write(commands1, file="download_MNase_files_1.sh")

quit()
n
chmod 770 download_MNase_files_1.sh
./download_MNase_files_1.sh

# cd /, R3.6.0
cell_line <- "HEPG2"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase/MNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access$score <- max(access$score) - access$score + min(access$score)
  access_data<-reduce(access[which(access$score >= quantile(access$score,(i)))])
  save(access_data,file=access_cell_line)
}




####### NOME #######
# NOMe-seq: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2387170
cd /storage/st20d/HiC_Kc167/alessandra/HEPG2/NOMe
fastq-dump SRR5008339
R3.5.0
cell_line <- "HEPG2" 
srr1 <- "SRR5008339"

commands1 <- c("#!/bin/bash\n\n")
cat <- paste0("cat ",srr1,".fastq > NOMe_",cell_line,".fastq", collapse =  "\n")
fastqc<-paste0("fastqc NOMe_",cell_line,".fastq -o fastqc_raw/")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U NOMe_",cell_line,".fastq -S 'NOMe_",cell_line,"_chip.sam' > 'NOMe_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t NOMe_",cell_line,"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n NOMe_",cell_line,"_peaks -B --nomodel > 'NOMe_",cell_line, "_peaksLog.log'")
commands1<-paste0(cat, "\n", fastqc, "\n", bowtie, "\n", macs2, "\n")
write(commands1, file="download_NOMe_files_1.sh")

quit()
n
chmod 770 download_NOMe_files_1.sh
./download_NOMe_files_1.sh

# cd /, R3.6.0
cell_line <- "HEPG2"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/NOMe"))
library(rtracklayer) 
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){ ### DOUBLE CHECKK THIS IS THE CORRECT FOR LOOP SINCE FOR IMR9P NOMe-seq I USED THE LIFTOVER.
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/NOMe/NOMe_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  save(access_data,file=access_cell_line)
}