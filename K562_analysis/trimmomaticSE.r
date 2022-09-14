rm(list=ls())

##get file names
filenames<-read.csv('/home/romana/Documents/HepG2/filenames_HepG2_ChIP_SE.txt',
          sep ='\t', head = F) ##for local

filenames$V1 <- gsub(".fastq.gz", "", filenames$V1)

#wtrite shell scripe
bash <- "#!/bin/bash\n"

for(i in filenames){
  bash <- c(bash, paste0("qsub ./trimmommaticWrapperSE.sh ", paste0(i, " ",
            i, "_trim ", "19" )))
}

write(bash, 'trim_HepG2_ChIP_SE.sh')
