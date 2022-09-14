################################################################################
##always use
################################################################################

rm(list=ls())

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

metadata <- read.csv('/home/romanapop/Dropbox/Documents/K562/SE/metadata_K562_ChIP_SE.tsv',
          sep ='\t', head = T) ##for local

##ctrl_meta <- read.csv('/home/romana/Documents/metadata_SE_ctrls.tsv',
          		##sep ='\t', head = T)  ##for local

##get file accession, TF name, treat, and ctrl columns
metadata <- metadata[,c(1,10,19, 47)]

##remove treated samples
metadata <- metadata[!grepl('a', as.character(metadata$Biosample.treatment),
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

##select rows with multiple controls
multiple <- metadata[grepl(',', metadata$Controlled.by),]

##separate the two files into different columns
multiple1 <- data.frame(do.call('rbind', strsplit(as.character(multiple$Controlled.by),
            ',',fixed=TRUE)))


multiple <- cbind(multiple, multiple1)

multiple$Controlled.by <- NULL

multiple.split <- split(multiple, as.character(multiple$Experiment.target))

##remove rows with multiple controls
metadata <- metadata[!grepl(',', metadata$Controlled.by),]

##remove rows with no control
metadata <- metadata[grepl('E', metadata$Controlled.by),]

##Create bash script to concatenate control files used for each exp
bash_commands <- c("#!/bin/bash\n",
                  "cd /home/rtpop/ENCODE_ChIP_data/HepG2/single_end/controls" )

for(i in seq_along(TF.split)){
    buffer <- TF.split[[i]][!duplicated(TF.split[[i]][,4]),] ##remove duplicate ctrl rows
    bash_commands<- c(bash_commands,paste0("cat ", paste0(buffer[,4],".fastq",
					collapse=" "),">",names(TF.split)[i],"_K562_SE_ctrl.fastq", collapse=""))
}

##add the exp that used multiple ctrls
for(i in seq_along(multiple.split)){
    buffer2 <- multiple.split[[i]][!duplicated(multiple.split[[i]][,4]),] ##remove duplicate ctrl rows
    bash_commands<- c(bash_commands,paste0("cat ", paste0(buffer2[,4],
      ".fastq ", buffer2[,5], ".fastq ", buffer2[,6], ".fastq ", buffer2[,7],
      ".fastq ", buffer2[,8], ".fastq ", buffer2[,9], ".fastq ", collapse=" "),
      ">",names(multiple.split)[i],
      "_K562_SE_ctrl.fastq", collapse=""))
}


write(bash_commands, 'cat_HepG2_ChIP_ctrl_SE.sh')

################################################################################


################################################################################
##only use when cat-ing experiment files
################################################################################

##Create bash script to concatenate files with data from the same TF
bash_commands <- c("#!/bin/bash\n",
                  "cd /home/rtpop/ENCODE_ChIP_data/K562/single_end/"
)

for(i in seq_along(TF.split)){
		bash_commands<- c(bash_commands,paste0("cat ", paste0(TF.split[[i]][,1],
                    ".fastq.gz", collapse=" "),">",names(TF.split)[i],
                    "_K562_SE.fastq.gz", collapse=""))
}

write(bash_commands, 'cat_K562_ChIP_SE.sh')

################################################################################
