rm(list=ls())

##get file names
filenames<-read.csv('/home/romana/Documents/K562/sam_filenames_K562_ChIP_ctrl_SE.txt',
          sep ='\t', head = F)

filenames$V1 <- gsub(".sam", "", filenames$V1)

#wtrite shell script
bash <- "#!/bin/bash\n"

for(i in filenames){
  bash <- c(bash, paste0("qsub ./sam2bamSEWrapper.sh ", paste0(i)))
}

write(bash, 'sam_to_bam_K562_ChIP_ctrl_SE.sh')
