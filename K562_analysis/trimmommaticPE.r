rm(list=ls())

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

##get file names
filenames<-read.csv('/home/romanapop/Dropbox/Documents/HepG2/filenames_HepG2_ChIP_PE.txt',
          sep ='_', head = F) ##for local

filenames$V5 <- gsub(".fastq.gz", "", filenames$V5)

##split by pair end
pair.split <- split(filenames, filenames$V1)

#wtrite shell scripe
bash <- "#!/bin/bash\n"

for(i in seq_along(pair.split)){
  bash <- c(bash, paste0("qsub ./trimmommaticWrapperPE.sh ", paste0(
          names(pair.split[i]), "_HepG2_ChIP_PE_1 ",
          names(pair.split[i]), "_HepG2_ChIP_PE_2 ",
          names(pair.split[i]), "_HepG2_ChIP_trim_PE_1 ",
          names(pair.split[i]), "_HepG2_ChIP_trim_unpaired_1 ",
          names(pair.split[i]), "_HepG2_ChIP_trim_PE_2 ",
          names(pair.split[i]), "_HepG2_ChIP_trim_unpaired_2 ", "19")))
}

write(bash, 'trim_HepG2_ChIP_PE.sh')
