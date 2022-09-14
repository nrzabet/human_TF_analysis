################### Pre-processing ATAC-seq data - general ####################

## Go to ENCODE https://www.encodeproject.org/search/?type=Experiment&status=released&assay_slims=DNA+accessibility&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&biosample_ontology.term_name=IMR-90&award.project=ENCODE&files.file_type=fastq&assay_title=ATAC-seq
## Select: DNA accessibility, ATAC-seq, Homo sapiens, IMR90 (or whatever cell line), ENCODE project, fastq files only, single ended if available
## Hit "download". Change "files.txt" to "files_ATAC.txt".
## Then open this file and copy-paste the first line into Google. It will download another file.
## Change "metadata.tsv" to "metadata_ATAC.tsv". 
## The metadata will be very short and it's normal: it's only one experiment per cell line, two replicates.
## Open cluster, screen, qrsh, go to IMR90 folder
mkdir -p ATAC
cd ATAC 
## Move files_ATAC.txt and metadata_ATAC.tsv to ATAC folder. 
## Go inside metadata_ATAC.tsv and change "IMR90" instead of "IMR-90" under "biosample term name". This will make things easier for future steps.
xargs -L 1 curl -O -J -L < files_ATAC.txt 
mkdir -p fastqc_raw

R3.6.0

metadata <- read.table("metadata_ATAC.tsv", sep="\t", header=T)

commands1 <- c("#!/bin/bash\n\n")
file <- metadata$"File.download.URL"
download <- paste0("wget ", file ,collapse =  "\n")
commands1 <- paste0(commands1, download,"\n")
commands1 <- paste0(commands1,"\n", "gunzip *.gz")
write(commands1, file="download_ATAC_files_1.sh")

commands2 <- c("#!/bin/bash\n\n")
splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")
commands2 <- paste0("cat ",splitfiles," > ATAC_",names(splitfiles),".fastq", collapse =  "\n")
write(commands2, file="download_ATAC_files_2.sh")

commands3 <- c("#!/bin/bash\n\n")
splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")
fastqc<-paste0("fastqc ATAC_",names(splitfiles),".fastq -o fastqc_raw/")
#after fastqc, check reads quality via the HTML page and decide whether to use cutadapt to remove some bases (based on "per base sequence content" module) and then run fastqc again, or move on to bowtie directly.
#in case I need to run cutadapt and fastqc2 here's the code (just need to change the number of bases to remove with -u argument)
cutadapt<-paste0("cutadapt -u 15 -o ATAC_", names(splitfiles),"_forward.fq.gz" , " ATAC_", names(splitfiles), ".fastq") # -u argument is the number of bases to delete. If it's positive (e.g. -u 5), it's the first 5 bases that are deleted. Otherwise e.g. (-u -5) it's the last 5 bases that are deleted. 
fastqc2<-paste0("fastqc ATAC_",names(splitfiles),"_forward.fq.gz")
bowtie<-paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U ATAC_",names(splitfiles),"_forward.fq.gz -S 'ATAC_",names(splitfiles),"_chip.sam' > 'ATAC_", names(splitfiles), "_bowtieLog.log'")
macs2<-paste0("macs2 callpeak -t ATAC_",names(splitfiles),"_chip.sam --broad --broad-cutoff 0.1 -g 2.7e9 -n ATAC_",names(splitfiles),"_peaks -B --nomodel > 'ATAC_", names(splitfiles), "_peaksLog.log'")
commands3<-paste0(fastqc, "\n", cutadapt, "\n", fastqc2, "\n", bowtie, "\n", macs2, "\n") # If I need to run cutadapt and fastqc2, use this line as "commands3"
#commands3<-paste0(fastqc, "\n", bowtie, "\n", macs2, "\n") # If I DON'T need to run cutadapt and fastqc2, use this line as "commands3"
write(commands3, file="download_ATAC_files_3.sh")

## IMPORTANT: ATACseq doesn't have controls, so for macs2 peakCalling I am not using controls like I did in the ChIP-seq data pre-processing.
## So -c argument has been deleted from macs2 callpeaks function.

quit()
n

chmod 770 download_ATAC_files_1.sh
chmod 770 download_ATAC_files_2.sh
chmod 770 download_ATAC_files_3.sh
./download_ATAC_files_1.sh
./download_ATAC_files_2.sh
./download_ATAC_files_3.sh # % overall alignment rate from BowtieLog 


### Now Generate QDA_IDs for ATACseq data to use for ChIPanalyser
# qrsh, cd /, R3.6.0
cell_line <- "IMR90"
setwd(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC"))
library(rtracklayer)
quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)
for(i in quantVec){
  access_cell_line <-paste0("access_",cell_line,"_",i,".Rda")
  access <- import(paste0("/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/ATAC/ATAC_",cell_line,"_peaks_treat_pileup.bdg"), format="bedGraph")
  access_data<-reduce(access[which(access$score >= quantile(access$score,i))])
  save(access_data,file=access_cell_line)
}

########################## PRE-PROCESSING DONE! ########################## Now move on to ChIPanalyser.