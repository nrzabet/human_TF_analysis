#####################################################
## Pre-processing of cell line, single ends #########
#####################################################
# This script has been generalised apart from Bash lines for which I've used IMR90 as example.

## Select your cell line from ENCODE, ChIP-seq only under "Assay title".
## Select fastq only from drop-down-menu of "Available file types", as well as single/paired ends based on needs.
## Hit download. Open "files.txt", run first line on Google. It will download a file called "metadata".
## Change metadata name to match cell line e.g. metadata_cell_line.tsv
## Do the exact same for ctrl files by selecting Controls ChIP-seq under "Assay title".
## Change files.txt to files_ctrl.txt, and metadata to metadata_ctrl_cell_line.tsv.
## So now you have 4 files: "files.txt", "files_ctrl.txt", "metadata_cell_line.tsv" and "metadata_ctrl_cell_line.tsv".
## Move everything to cluster. Open screen, qrsh, go to right directory (HiC_Kc167/alessandra/cell_line/preProcessing/)
xargs -L 1 curl -O -J -L < files.txt 
xargs -L 1 curl -O -J -L < files_ctrl.txt # It'll download fastq files (will take a while)
## Open "metadata_cell_line.tsv" and go all the way to "Controlled.by" column.
## There will only be a number of files, in this case it's only 4. 
## So now from all the downloaded fastq control files (line 15 of this script), I can manually remove the ones I do not need.

## In the same working directory (alessandra/cell_line/preProcessing/), type:
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
## Now go inside Trimmomatic-0.39 and move the whole content outside of its folder, into preProcessing folder.

########################## Concatenating ##########################
R3.5.0
## ChIP data
rm(list=ls())
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

cell_line <- "IMR90" # change according to needs

metadata <-read.csv(paste0("metadata_",cell_line,".tsv"),sep ='\t', head = T) # Get file names for ChIP data
metadata <- metadata[,c(1,10,22,50)] # Get file accession, TF name, treatments, and ctrl columns
metadata$Controlled.by <- gsub("files", "", metadata$Controlled.by)
metadata$Controlled.by <- gsub('/', '', metadata$Controlled.by) # Remove extra stuff from file name
# Now by typing metadata$Controlled.by, I can see I have 4 types of control files:
# in this case they are ENCFF000YCQ, ENCFF000YCS, ENCFF584LXI and ENCFF492DKQ.
TF.split <- split(metadata, metadata$Experiment.target) # Split by TF

# Create bash script to concatenate files with data from the same TF
bash_commands <- c("#!/bin/bash\n",
                   paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in seq_along(TF.split)){
  bash_commands<- c(bash_commands,paste0("cat ", paste0(TF.split[[i]][,1],
                                                        ".fastq.gz", collapse=" "),
                                         ">",names(TF.split)[i],"_",cell_line,"_SE.fastq.gz", collapse=""))
}
write(bash_commands, paste0("cat_",cell_line,"_ChIP_SE.sh"))

## Controls
rm(list=ls())

metadata<-read.table(paste0("metadata_",cell_line,".tsv"), sep="\t", header=T, stringsAsFactors = FALSE) # Same metadata as ChIP file, the one from Encode.
TF_filenames <- unlist(lapply(lapply(split(metadata$File.accession,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),paste0,".fastq"), paste0, collapse=" ")) # Extract TF ChIP files
control_filenames_clean <- as.character(metadata$Controlled.by)
control_filenames_clean<-gsub(",", "", control_filenames_clean)
control_filenames_clean<-gsub("/files/", "", control_filenames_clean)
control_filenames_clean<-gsub("/", "_", control_filenames_clean)
control_filenames_clean<-gsub(" ", "", control_filenames_clean)
control_filenames <- unlist(lapply(lapply(split(control_filenames_clean,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),unique), paste0, collapse="_"))
control_filenames <- gsub("__", "_", control_filenames)
control_filenames_ordered <- strsplit(control_filenames, "_")
control_filenames_ordered <- lapply(control_filenames_ordered, sort)
control_filenames_ordered <- lapply(control_filenames_ordered, paste0, "_")
control_filenames <- sapply(control_filenames_ordered, paste0, collapse="")
control_filenames_unique <- unique(control_filenames) # Extract control ChIP files

commands_controls <- c("#!/bin/bash\n\n")
cat <- paste0("cat ",gsub("_",".fastq.gz ", control_filenames_unique)," > ChIPseq_control_",control_filenames_unique,"raw.fastq.gz")
trim<- paste0("java -jar trimmomatic-0.39.jar SE -phred33 ChIPseq_control_",control_filenames_unique,"raw.fastq.gz ChIPseq_control_",control_filenames_unique,"trimmed.fq.gz  ILLUMINACLIP:/storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25")
bowtie2<- paste0("bowtie2 -x /storage/st05d/ReferenceGenomes/Bowtie2Indexes/hg38 -U ChIPseq_control_",control_filenames_unique,"trimmed.fq.gz -S ChIPseq_control_",control_filenames_unique,"chip.sam")
commands_controls<- c(commands_controls, paste0(cat, "\n", trim, "\n", bowtie2, "\n"))
write(commands_controls, paste0(cell_line,"_control_files.sh"))

quit()

chmod 770 cat_IMR90_ChIP_SE.sh
chmod 770 IMR90_control_files.sh
./cat_IMR90_ChIP_SE.sh
./IMR90_control_files.sh # Executing .sh files

########################## Fastqc ##########################
## Create a .txt file (notepad) with all the fastq.gz file names (they look like: TFname-human_cellline_SE.fastq.gz).
## Save as filenames_cellline. Do it in the WinSC.
## Do the same for control files and save as filenames_cellline_ctrl. 
R3.5.0
# ChIP data:
rm(list=ls())
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

filenames <- read.csv(paste0("filenames_",cell_line,".txt"), sep ="", head = F)
bash_commands <- c("#!/bin/bash\n",
                   paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames$V1){
  bash_commands <- c(bash_commands, paste0("qsub ./fastqcWrapper.sh ",
                                           paste0(i, collapse = " ")))
}
write(bash_commands, paste0("fastqc_",cell_line,"_ChIP_SE.sh"))

# Control data:
rm(list=ls())
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

filenames <- read.csv(paste0("filenames_",cell_line,"_ctrl.txt"), sep ="", head = F)
bash_commands <- c("#!/bin/bash\n",
                   paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames$V1){
  bash_commands <- c(bash_commands, paste0("qsub ./fastqcWrapper.sh ",
                                           paste0(i, collapse = " ")))
}
write(bash_commands, paste0("fastqc_",cell_line,"_ctrl_SE.sh"))

quit()

## Move "fastqcWrapper.sh" to "alessandra/cell_line/preProcessing" folder in the cluster.
chmod 770 fastqc_IMR90_ChIP_SE.sh 
chmod 770 fastqc_IMR90_ctrl_SE.sh
chmod 770 fastqcWrapper.sh
exit
mkdir -p fastqc
./fastqc_IMR90_ChIP_SE.sh
./fastqc_IMR90_ctrl_SE.sh
## It'll say that jobs have been submitted. Type qstat -u ap16539 to know when it's done (if doesn't return anything).
## If it says "no such file or directory", it's common when creating a .txt file on Windows, which was then moved to Linux.
## So just re-create filenames_celline.txt or fastqc_....sh from WinSCP.

########################## Trimmomatic ##########################
## Screen, qrsh, redirect from pwd to preProcessing folder.
R3.5.0
# ChIP data
rm(list=ls())
filenames <- read.csv(paste0("filenames_",cell_line,".txt"), sep ='\t', head = F)
filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

bash <- c("#!/bin/bash\n",
          paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames){
  bash <- c(bash, paste0("qsub ./trimmommaticWrapperSE.sh ", paste0(i, " ",
                                                                    i, "_trim ", "19" )))
}
write(bash, paste0("trim_",cell_line,"_ChIP_SE.sh"))

# Control data
rm(list=ls())
filenames <- read.csv(paste0("filenames_",cell_line,"_ctrl.txt"), sep ='\t', head = F)
filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

bash <- c("#!/bin/bash\n",
          paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames){
  bash <- c(bash, paste0("qsub ./trimmommaticWrapperSE.sh ", paste0(i, " ",
                                                                    i, "_trim ", "19" )))
}
write(bash, paste0("trim_",cell_line,"_ctrl_SE.sh"))

quit()

## Move "trimmmaticWrapperSE.sh" to "alessandra/cell_line/preProcessing" folder in the cluster.
chmod 770 trim_IMR90_ChIP_SE.sh
chmod 770 trim_IMR90_ctrl_SE.sh
chmod 770 trimmommaticWrapperSE.sh
exit
./trim_IMR90_ChIP_SE.sh
./trim_IMR90_ctrl_SE.sh

########################## Bowtie ##########################
## Create a .txt file (notepad) with all the trim_fastq.gz file names (they look like: TFname-human_cellline_SE_trim.fastq.gz).
## Save as trim_filenames_cellline. Do it in the WinSC.
## Do the same for control files and save as trim_filenames_cellline_ctrl. 
# Qrsh, then redirect from pwd
R3.5.0

# ChIP data
rm(list=ls())
filenames<-read.csv(paste0("trim_filenames_",cell_line,".txt"), sep ='\t', head = F)
filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

bash <- c("#!/bin/bash\n",
          paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames){
  bash <- c(bash, paste0("qsub ./bowtieWrapperSE.sh ", paste0(i, " ", i, "_alignment_stats")))
}
write(bash, paste0("align_",cell_line,"_ChIP_SE.sh"))

# Control data
rm(list=ls())
filenames<-read.csv(paste0("trim_filenames_",cell_line,"_ctrl.txt"), sep ='\t', head = F)
filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

bash <- c("#!/bin/bash\n",
          paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in filenames){
  bash <- c(bash, paste0("qsub ./bowtieWrapperSE.sh ", paste0(i, " ", i, "_alignment_stats")))
}
write(bash, paste0("align_",cell_line,"_ctrl_SE.sh"))

quit()

## Move "bowtieWrapperSE.sh" to "alessandra/cell_line/preProcessing" folder in the cluster.
chmod 770 align_IMR90_ChIP_SE.sh 
chmod 770 align_IMR90_ctrl_SE.sh
chmod 770 bowtieWrapperSE.sh 
exit
./align_IMR90_ChIP_SE.sh
./align_IMR90_ctrl_SE.sh


########################## Peak calling ##########################
## Create a .txt file (notepad) with all the sam file names (they look like: TFname-human_cellline_SE_trim.sam).
## Save as sam_filenames_cellline. Do it in the WinSC.
## Do the same for control files and save as sam_filenames_cellline_ctrl.
R3.5.0
rm(list=ls())
options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

filenames <- read.csv(paste0("sam_filenames_",cell_line,".txt"), sep ='\t', head = F)
ctrl_filenames <- read.csv(paste0("sam_filenames_",cell_line,"_ctrl.txt"), sep ='\t', head = F)
filenames<-cbind(filenames, ctrl_filenames)
colnames(filenames) <- c("ChIP", "ctrls")
filenames$ChIP <- gsub(".sam", "", filenames$ChIP)
filenames$ctrls <- gsub(".sam", "", filenames$ctrls)

bash <- c("#!/bin/bash\n",
          paste0("cd /storage/st20d/HiC_Kc167/alessandra/",cell_line,"/preProcessing"))
for(i in rownames(filenames)){
  bash <- c(bash, paste0("qsub ./callpeaksWrapper.sh ", paste0(filenames[i,1],
                                                               " ", filenames[i,2])))
}
write(bash, paste0("callpeaks_",cell_line,"_ChIP_ctrl_SE.sh"))

quit()

## Move "callPeaksWrapperSE.sh" to "alessandra/cell_line/preProcessing" folder in the cluster.
chmod 770 callpeaks_IMR90_ChIP_ctrl_SE.sh 
chmod 770 callpeaksWrapper.sh 
exit 
./callpeaks_IMR90_ChIP_ctrl_SE.sh

########################## PRE-PROCESSING DONE! ##########################