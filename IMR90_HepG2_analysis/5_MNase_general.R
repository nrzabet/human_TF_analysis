################### Pre-processing MNase-seq data - general ####################

## Go to this page https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44985 and click on SRA number "SRA SRP019045" at the bottom.
## A new page will pop up. Go to the bottom to find the SRR numbers. Write these numbers down: SRR769551, SRR769552, SRR769553.
## Double check if these reads are single or paired end reads (for IMR90 these should be single).
## Open cluster, screen, qrsh, go to IMR90 folder
mkdir -p MNase
cd MNase
mkdir -p fastqc_raw
## Now run fastq-dump on the 3 SRR numbers above. This will automatically download the fastq files in the cluster.
fastq-dump SRR769551
fastq-dump SRR769552
fastq-dump SRR769553 # this is the line to run for single-ends. For paired-ends check "fastq-dump -I --split-files SRR769551"

R3.5.0
cell_line <- "IMR90" # change according to needs
srr1 <- "SRR769551"
srr2 <- "SRR769552"
srr3 <- "SRR769553" # change these numbers according to needs with the different SRR numbers. Add more if necessary.

commands1<- c("#!/bin/bash\n\n")
commands1 <- paste0("cat ",srr1,".fastq ",srr2,".fastq ",srr3,".fastq"," > MNase_",cell_line,".fastq", collapse =  "\n") # add more "srr1" etc according to needs.
write(commands1, file="download_MNase_files_1.sh") # this will concatenate the SRR fastq files into one

commands2 <- c("#!/bin/bash\n\n")
fastqc<-paste0("fastqc MNase_",cell_line,".fastq -o fastqc_raw/")
#after fastqc, check reads quality via the HTML page and decide whether to use cutadapt to remove some bases(based on "per base sequence content" module) and then run fastqc again, or move on to bowtie directly.
#in case I need to run cutadapt and fastqc2 here's the code (just need to change the number of bases to remove with -u argument)
cutadapt<-paste0("cutadapt -u 4 -o MNase_",cell_line,"_forward.fq.gz"," MNase_",cell_line,".fastq") # -u argument is the number of bases to delete. If it's positive (e.g. -u 5), it's the first 5 bases that are deleted. Otherwise e.g. (-u -5) it's the last 5 bases that are deleted. 
fastqc2<-paste0("fastqc MNase_",cell_line,"_forward.fq.gz")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U MNase_",cell_line,"_forward.fq.gz -S 'MNase_",cell_line,"_chip.sam' > 'MNase_",cell_line,"_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t MNase_",cell_line,"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n MNase_",cell_line,"_peaks -B --nomodel > 'MNase_",cell_line, "_peaksLog.log'")
commands2<-paste0(fastqc, "\n", cutadapt, "\n", fastqc2, "\n", bowtie, "\n", macs2, "\n") # If I need to run cutadapt and fastqc2, use this line as "commands2"
#commands2<-paste0(fastqc, "\n", bowtie, "\n", macs2, "\n") # If I DON'T need to run cutadapt and fastqc2, use this line as "commands2"
write(commands2, file="download_MNase_files_2.sh")

## IMPORTANT: MNase doesn't have controls, so for macs2 peakCalling I am not using controls like I did in the ChIP-seq data pre-processing.
## So -c argument has been deleted from macs2 callpeaks function.

quit()
n

chmod 770 download_MNase_files_1.sh
chmod 770 download_MNase_files_2.sh
./download_MNase_files_1.sh
./download_MNase_files_2.sh


### Now Generate QDA_IDs for MNase data to use for ChIPanalyser
# qrsh, cd /, R3.5.0
cell_line <- "IMR90"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/MNase/MNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access$score <- max(access$score) - access$score + min(access$score)
  access_data<-reduce(access[which(access$score >= quantile(access$score,(i)))]) ## THIS LINE IS IMPORTANT FOR MNase-seq.
  save(access_data,file=access_cell_line)
}

########################## PRE-PROCESSING DONE! ########################## Now move on to ChIPanalyser.