rm(list=ls())

##get file names
filenames <- read.csv('/home/romanapop/Dropbox/Documents/HepG2/filenames_HepG2_ChIP_PE.txt',
            sep ='', head = F) ##for local

##generate bash script
bash_commands <- "#!/bin/bash\n"

for(i in filenames$V1){
  bash_commands <- c(bash_commands, paste0("qsub ./fastqcWrapper.sh ", paste0(i,
                    collapse = " ")))
}

write(bash_commands, 'fastqc_HepG2_ChIP_PE.sh')
