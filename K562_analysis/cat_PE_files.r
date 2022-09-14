################################################################################
##always use
################################################################################
rm(list=ls())

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

##get file names
metadata<-read.csv('/home/romanapop/Dropbox/Documents/HepG2/metadata_HepG2_ChIP_ctrl_PE.tsv',
          sep ='\t', head = T) ##for local

##selecting columns of interest
metadata <- metadata[,c(1,10, 19, 36, 47)]

##remove treated files
metadata <- metadata[!grepl('inter', as.character(metadata$Biosample.treatment),
        ignore.case =T),]

##remove extra stuff from file name
metadata$Controlled.by <- gsub("files", "", metadata$Controlled.by)
metadata$Controlled.by <- gsub('/', '', metadata$Controlled.by)

##split by TF
TF.split <- split(metadata, as.character(metadata$Experiment.target))
################################################################################

################################################################################
##only use when cat-ing ctrl files
################################################################################
multiple <- metadata[grepl(',', as.character(metadata$Controlled.by)),]

##separate the two files into different columns
multiple1 <- data.frame(do.call('rbind', strsplit(as.character(multiple$Controlled.by),
            ',',fixed=TRUE)))

multiple <- cbind(multiple, multiple1)

multiple$Controlled.by <- NULL

multiple.split <- split(multiple, as.character(multiple$Experiment.target))

##remove rows with multiple controls
metadata <- metadata[!grepl(',', as.character(metadata$Controlled.by)),]

##remove rows with no control
metadata <- metadata[grepl('E', as.character(metadata$Controlled.by)),]

##split by TF after removing multiple ctrls and no ctrls
TF.split <- split(metadata, as.character(metadata$Experiment.target))

##Create bash script to concatenate control files used for each exp
bash_commands <- "#!/bin/bash\n"
bash_commands <- c(bash_commands, "cd /home/rtpop/ENCODE_ChIP_data/K562/pair_end/controls\n",
		 "gunzip *.fastq.gz\n")

for(i in seq_along(TF.split)){
   buffer <- TF.split[[i]][!duplicated(TF.split[[i]][,c(4,5)]),] ##remove duplicate ctrl rows
   buffer.split <- split(buffer, buffer$Paired.end)
      for(j in seq_along(buffer.split)){
        if(j == 1){
         bash_commands<- c(bash_commands,paste0("cat ",paste0(buffer.split[[j]][,5][!duplicated(buffer.split[[j]][,5])],
           ".fastq", collapse=" "),">",
           names(TF.split)[i], "_HepG2_ChIP_ctrl_PE_1.fastq", collapse=""))
    } else{
        bash_commands<- c(bash_commands,paste0("cat ", paste0(buffer.split[[j]][,5][!duplicated(buffer.split[[j]][,5])],
          ".fastq", collapse=" "),">",names(TF.split)[i],
          "_HepG2_ChIP_ctrl_PE_2.fastq", collapse=""))
    }
  }
}

##add the exp that used multiple ctrls
for(i in seq_along(multiple.split)){
    buffer3 <- multiple.split[[i]][!duplicated(multiple.split[[i]][,c(4,5,6)]),] ##remove duplicate ctrl rows
      for(j in buffer3[,4]){
        if(j == 1){
          bash_commands<- c(bash_commands,paste0("cat ",paste0(as.character(buffer3[1,5]),
            ".fastq", as.character(buffer3[1,6]), collapse=" "),">",
            names(multiple.split)[i], "_HepG2_ChIP_ctrl_PE_1.fastq", collapse=""))
        } else{
          bash_commands<- c(bash_commands,paste0("cat ", paste0(as.character(buffer3[1,5]),
            ".fastq", as.character(buffer3[1,6]), collapse=" "),">",names(multiple.split)[i],
            "_HepG2_ChIP_ctrl_PE_2.fastq", collapse=""))
        }
    }
}

bash_commands <- c(bash_commands, "gzip *_HepG2_SE_ctrl.fastq\n")

write(bash_commands, 'cat_HepG2_ChIP_ctrls_PE.sh')

################################################################################


################################################################################
##only use when cat-ing experiment files
################################################################################

##Create bash script to concatenate files
bash_commands <- "#!/bin/bash\n"
bash_commands <- c(bash_commands, "cd /home/rtpop/ENCODE_ChIP_data/HepG2/pair_end/controls\n",
                  "gunzip *fastq.gz")

for(i in seq_along(TF.split)){
   buffer <- TF.split[[i]]
   buffer.split <- split(buffer, buffer$Paired.end)
      for(j in seq_along(buffer.split)){
        if(j == 1){
         bash_commands<- c(bash_commands,paste0("cat ",paste0(as.character(buffer.split[[j]][,1]),
           ".fastq", collapse=" "),">",
           names(TF.split)[i], "_HepG2_ChIP_PE_1.fastq", collapse=""))
    } else{
        bash_commands<- c(bash_commands,paste0("cat ", paste0(as.character(buffer.split[[j]][,1]),
          ".fastq", collapse=" "),">",names(TF.split)[i],
          "_HepG2_ChIP_PE_2.fastq", collapse=""))
    }
  }
}

bash_commands <- c(bash_commands, "gzip *.fastq\n")

write(bash_commands, 'cat_HepG2_ChIP_PE.sh')
################################################################################
