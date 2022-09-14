################### Pre-processing DNA-seq data (DNase I) - general ####################

## Go to ENCODE https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=DNase-seq&assay_slims=DNA+accessibility&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=IMR-90&award.project=ENCODE&files.file_type=fastq
## Select: DNA accessibility, DNase-seq, Homo sapiens, IMR90 (or whatever cell line), ENCODE project, fastq files only, single ended
## Hit "download". Change "files.txt" to "files_DNase.txt".
## Then open this file and copy-paste the first line into Google. It will download another file.
## Change "metadata.tsv" to "metadata_DNase.tsv". 
## The metadata will be very short and it's normal: it's only one experiment per cell line, two replicates.
## Open cluster, screen, qrsh, go to IMR90 folder
mkdir -p DNase
cd DNase 
## Move files_DNase.txt and metadata_DNase.tsv to DNase folder. 
## Go inside metadata_DNase.tsv and change "IMR90" instead of "IMR-90" under "biosample term name". This will make things easier for future steps.
xargs -L 1 curl -O -J -L < files_DNase.txt 
mkdir -p fastqc_raw

R3.5.0

metadata <- read.table("metadata_DNase.tsv", sep="\t", header=T)

commands1 <- c("#!/bin/bash\n\n")
file <- metadata$"File.download.URL"
download <- paste0("wget ", file ,collapse =  "\n")
commands1 <- paste0(commands1, download,"\n")
commands1 <- paste0(commands1,"\n", "gunzip *.gz")
write(commands1, file="download_DNase_files_1.sh")

commands2 <- c("#!/bin/bash\n\n")
splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")
commands2 <- paste0("cat ",splitfiles," > DNase_",names(splitfiles),".fastq", collapse =  "\n")
write(commands2, file="download_DNase_files_2.sh")

commands3 <- c("#!/bin/bash\n\n")
splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")
fastqc<-paste0("fastqc DNase_",names(splitfiles),".fastq -o fastqc_raw/")
cutadapt<-paste0("cutadapt -u -30 -o DNase_", names(splitfiles),"_forward.fq.gz" , " DNase_", names(splitfiles), ".fastq") # -u argument is the number of bases to delete. If it's positive (e.g. -u 5), it's the first 5 bases that are deleted. Otherwise e.g. (-u -5) it's the last 5 bases that are deleted. 
#after fastqc, check reads quality via the HTML page and decide whether to use cutadapt to remove some bases (based on "per base sequence content" module) and then run fastqc again, or move on to bowtie directly.
#in case I need to run cutadapt and fastqc2 here's the code (just need to change the number of bases to remove with -u argument)
#trim<- paste0("java -jar trimmomatic-0.39.jar SE -threads 5 DNase_",names(splitfiles),".fastq DNase_",names(splitfiles),"_forward.fq.gz 
# ILLUMINACLIP:/storage/st20d/HiC_Kc167/alessandra/IMR90/preProcessing/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25") # In this case trimmommatic was deleted, instead I used cutadapt
fastqc2<-paste0("fastqc DNase_",names(splitfiles),"_forward.fq.gz")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U DNase_",names(splitfiles),"_forward.fq.gz -S 'DNase_",names(splitfiles),"_chip.sam' > 'DNase_", names(splitfiles), "_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t DNase_",names(splitfiles),"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n DNase_",names(splitfiles),"_peaks -B --nomodel > 'DNase_", names(splitfiles), "_peaksLog.log'")
commands3<-paste0(fastqc, "\n", cutadapt, "\n", fastqc2, "\n", bowtie, "\n", macs2, "\n") # If I need to run cutadapt and fastqc2, use this line as "commands3"
#commands3<-paste0(fastqc, "\n", bowtie, "\n", macs2, "\n") # If I DON'T need to run cutadapt and fastqc2, use this line as "commands3"
write(commands3, file="download_DNase_files_3.sh")

## IMPORTANT: DNase doesn't have controls, so for macs2 peakCalling I am not using controls like I did in the ChIP-seq data pre-processing.
## So -c argument has been deleted from macs2 callpeaks function.

quit()
n

chmod 770 download_DNase_files_1.sh
chmod 770 download_DNase_files_2.sh
chmod 770 download_DNase_files_3.sh
./download_DNase_files_1.sh
./download_DNase_files_2.sh
./download_DNase_files_3.sh # 98.01% overall alignment rate from BowtieLog 


### Now Generate QDA_IDs for DNase data to use for ChIPanalyser
# qrsh, cd /, R3.5.0
cell_line <- "IMR90"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/DNase/DNase_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  save(access_data,file=access_cell_line)
}

########################## PRE-PROCESSING DONE! ########################## Now move on to ChIPanalyser.