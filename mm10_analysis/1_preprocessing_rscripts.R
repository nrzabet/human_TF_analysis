##########
#DNase

	setwd("/storage/projects/ZabetLab/mm10_ENCODE/data/")

	metadata<-read.table("../scripts/metadata_DNase.tsv", sep="\t", header=T)

	#download	

	commands1 <- c("#!/bin/bash\n\n")
		
		file<- metadata$"File.download.URL"

		download<-paste0("wget ", file ,collapse =  "\n")


	commands1<- paste0(commands1, download,"\n")
	commands1<- paste0(commands1,"\n", "gunzip *.gz")

	write(commands1, file="../scripts/download_DNase_files_1.sh")

	#cat

	commands2 <- c("#!/bin/bash\n\n")

		splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")

	commands2 <- paste0("cat ",splitfiles," > DNase_",names(splitfiles),".fastq", collapse =  "\n")
	write(commands2, file="../scripts/download_DNase_files_2.sh")

	#preprocess

	commands3 <- c("#!/bin/bash\n\n")
		splitfiles<- lapply(lapply(split(metadata$File.accession, metadata$Biosample.term.name), paste0, ".fastq"), paste0, collapse=" ")

		
		fastqc<- paste0("fastqc DNase_",names(splitfiles),".fastq -o fastqc_raw/")
		trim<- paste0("java -jar trimmomatic-0.39.jar SE -threads 5 DNase_",names(splitfiles),".fastq DNase_",names(splitfiles),"_forward.fq.gz  ILLUMINACLIP:/home/dn19121/project/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25")
		fastqc2<-paste0("fastqc DNase_",names(splitfiles),"_forward.fq.gz")
		bowtie<- paste0("bowtie2 -x /storage/projects/ReferenceGenomes/Bowtie2Indexes/mm10 -U DNase_",names(splitfiles),"_forward.fq.gz -S DNase_",names(splitfiles),"_chip.sam")
		macs2<-paste0("macs2 callpeak -t DNase_",names(splitfiles),"_chip.sam --broad --broad-cutoff 0.1 -g 1.87e9 -n DNase_",names(splitfiles),"_peaks -B")


	commands3<- paste0(fastqc, "\n", trim, "\n", fastqc2, "\n", bowtie, "\n", macs2, "\n")
	write(commands3, file="../scripts/download_DNase_files_3.sh")


###########
#control

	metadata<-read.table("metadata.tsv", sep="\t", header=T, stringsAsFactors = FALSE)

	#extract TF ChIP files
	TF_filenames <- unlist(lapply(lapply(split(metadata$File.accession,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),paste0,".fastq"), paste0, collapse=" "))


	#extract control ChIP files
	control_filenames_clean <- as.character(metadata$Controlled.by)
	control_filenames_clean<-gsub(",", "", control_filenames_clean)
	control_filenames_clean<-gsub("/files/", "", control_filenames_clean)
	control_filenames_clean<-gsub("/", "_", control_filenames_clean)
	control_filenames_clean<-gsub(" ", "", control_filenames_clean)

	control_filenames <- unlist(lapply(lapply(split(control_filenames_clean,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),unique), paste0, collapse="_"))
	control_filenames <- gsub("__", "_", control_filenames)
	control_filenames_unique <- unique(control_filenames)

	#preprocess_controls

	commands_controls <- c("#!/bin/bash\n\n")

		cat <- paste0("cat ",gsub("_",".fastq ", control_filenames_unique)," > ChIPseq_control_",control_filenames_unique,"raw.fastq")
		trim<- paste0("java -jar trimmomatic-0.39.jar SE -threads 5 ChIPseq_control_",control_filenames_unique,"raw.fastq ChIPseq_control_",control_filenames_unique,"trimmed.fq.gz  ILLUMINACLIP:/home/dn19121/project/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25")
		bowtie2<- paste0("bowtie2 -x /storage/projects/ReferenceGenomes/Bowtie2Indexes/mm10 -U ChIPseq_control_",control_filenames_unique,"trimmed.fq.gz -S ChIPseq_control_",control_filenames_unique,"chip.sam")

	commands_controls<- c(commands_controls, paste0(cat, "\n", trim, "\n", bowtie, "\n"))
	write(commands_controls, "preprocess_control_files.sh")


############
#TF data

	setwd("/storage/projects/ZabetLab/mm10_ENCODE/data/")

	metadata<-read.table("../scripts/metadata.tsv", sep="\t", header=T)

	commands1 <- c("#!/bin/bash\n\n")

		for(i in 1:nrow(metadata)){

			download<-list()
			
			#download with correct output name
				download[as.character(i)]<- paste0("wget ", metadata[i,"File.download.URL"])

			download<-unlist(download)
			
			commands1<- paste0( download,"\n")
		}

		commands1<- paste0(commands1,"\n", "gunzip *.gz")
	write(commands1, "../scripts/download_TF_files_1.sh")

	commands2 <- c("#!/bin/bash\n\n")

		split_filenames <- unlist(lapply(lapply(split(metadata$File.accession,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),paste0,".fastq"), paste0, collapse=" "))

		commands2 <- paste0("cat ",split_filenames," > ChIPseq_",names(split_filenames),".fastq")
	write(commands2, "../scripts/download_TF_files_2.sh")

	commands3 <- c("#!/bin/bash\n\n")

		split_filenames <- unlist(lapply(lapply(split(metadata$File.accession,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),paste0,".fastq"), paste0, collapse=" "))
		#names(split_filenames)
		
		
		fastqc<- paste0("fastqc ChIPseq_",names(split_filenames),".fastq -o fastqc_raw/")
		trim<- paste0("java -jar trimmomatic-0.39.jar SE -threads 5 ChIPseq_",names(split_filenames),".fastq ChIPseq_",names(split_filenames),"_forward.fq.gz  ILLUMINACLIP:/home/dn19121/project/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32")
		fastqc2<-paste0("fastqc ChIPseq_",names(split_filenames),"_forward.fq.gz")
		bowtie<- paste0("bowtie2 -x /storage/projects/ReferenceGenomes/Bowtie2Indexes/mm10 -U ChIPseq_",names(split_filenames),"_forward.fq.gz -S ChIPseq_",names(split_filenames),"_chip.sam")
		
	commands3<- paste0(fastqc, "\n",trim, "\n", fastqc2, "\n",bowtie, "\n")
	write(commands3, "../scripts/download_TF_files_3.sh")


	commands_macs2 <- c("#!/bin/bash\n\n")
		TF_filenames <- unlist(lapply(lapply(split(metadata$File.accession,paste0(metadata$Biosample.term.name,"_",metadata$Experiment.target)),paste0,".fastq"), paste0, collapse=" "))

		macs2<-paste0("macs2 callpeak -t ChIPseq_",names(TF_filenames),"_chip.sam -c ChIPseq_control_", control_filenames,"chip.sam -f SAM -g 1.87e9 -n ChIPseq_",names(TF_filenames),"_peaks -B -q 0.05")

	commands_macs2<- c(commands_macs2, paste0(macs2, "\n"))
	write(commands_macs2, "preprocess_call_peaks.sh")





#executing: set the directory 
#chmod +x ../../scripts/script-name-here.sh
#bash ../../scripts/script-name-here.sh


