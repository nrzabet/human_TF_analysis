rm(list=ls())

##get filenames of sam files
filenames <- read.csv('/home/romanapop/Dropbox/Documents/HepG2/K562_sam_filenames_ChIP_CTRL_SE.txt',
          sep ='\t', head = F)

ctrl_filenames <- read.csv('/home/romana/Documents/K562/sam_filenames_K562_ChIP_ctrl_PE_nomodel.txt',
          #sep ='\t', head = F)

##get ones that had no controls specified
no_ctrl<-filenames[which(!filenames$V1 %in% ctrl_filenames$V1),1]
filenames<-cbind(filenames, ctrl_filenames)
colnames(filenames) <- c("ChIP", "ctrls")
filenames$ChIP <- gsub(".sam", "", filenames$ChIP)
filenames$ctrls <- gsub(".sam", "", filenames$ctrls)

#wtrite shell script
bash <- "#!/bin/bash\n"

for(i in rownames(filenames)){
  bash <- c(bash, paste0("qsub ./callpeaksWrapper.sh ", paste0(filenames[i,1],
          " ", filenames[i,2])))
}

write(bash, 'callpeaks_K562_ChIP_ctrl_SE.sh')
