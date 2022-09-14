rm(list=ls())

#get filenames of trimmed fastq files
filenames<-read.csv('/home/romana/Documents/K562/trim_filenames_K562_ChIP_SE.txt',
          sep ='\t', head = F) ##for local
filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

#wtrite shell script
bash <- "#!/bin/bash\n"

for(i in filenames){
  bash <- c(bash, paste0("qsub ./bowtieWrapperSE.sh ", paste0(i, " ", i, "_alignment_stats")))
}

write(bash, "align_K562_ChIP_SE.sh")
