rm(list=ls())

#getting the filenames for trimmed fastq files
filenames <- read.csv('/home/romanapop/Dropbox/Documents/K562/trimmed_filenames_K562_ChIP_PE.txt',
          sep ='_', head = F) ##for local
filenames$V5 <- gsub(".fastq.gz", "", filenames$V5)

#splitting by TF
TF.split <- split(filenames, filenames$V1)

#wtrite shell script
bash <- "#!/bin/bash\n"

for(i in seq_along(TF.split)){
  bash <- c(bash, paste0("qsub ./bowtieWrapperPE.sh ", paste0(names(TF.split[i]),
            "_K562_ChIP_trim_PE_1 ", paste0(names(TF.split[i]),
            "_K562_ChIP_trim_PE_2 ", names(TF.split[i]),
            "_K562_PE"))))
}

write(bash, "align_K562_ChIP_PE.sh")
