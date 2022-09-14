#table for comparing cell lines
#matrices:cell lines
#inner matrices: transcription factors
#rows:quantvec
#cols bound molecules lamdba optimal auc

setwd("/storage/projects/ZabetLab/mm10_ENCODE/scripts/")
library(GenomicRanges)
library(ChIPanalyser)

quantVec<-c(seq(0.0,0.9, by=0.1), 0.95,0.99)

############################
#function to get data into a dataframe
#with corresponding TF names

AUC_optimal <- rep(0, length(quantVec))
lambda <- rep(0, length(quantVec))
bm <- rep(0, length(quantVec))

getCellLineData<- function (TF){
	for(QDA_ID in 1:length(quantVec)){
		if (file.exists(paste0("../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))){
		load(paste0("../objects/set2/optimal_validation_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
		AUC_optimal[QDA_ID] <- optimal_validation_TF$GOF@profiles[[1]][[1]]["AUCMean"]
		
		load(paste0("../objects/set2/optimal_",cell_line,"_",TF,"_",quantVec[QDA_ID],"_gw.Rda"))
		param <- optimal_TF[[1]][[1]][["MSE"]]
		lambda[QDA_ID]<- param[1]
		bm[QDA_ID] <- param[2]
		}else{
		param<- "NA"
		lambda[QDA_ID]<- "NA"
		bm[QDA_ID] <- "NA"
		AUC_optimal[QDA_ID]<-"NA"
		}
		
		}
	df <- data.frame( AUC_optimal,lambda, bm, stringsAsFactors = F)
	rownames(df)<- c(quantVec)
	colnames(df)<- c( "AUC", "lambda", "bound mol")
	#how to name the elements of the list as TF names?
	#names(df)<- c(TF)
	return(df)
	}
	
############################
###analyse
metadata<-read.table("../scripts/metadata.tsv", sep="\t", header=T)

split<-split(metadata$Experiment.target,metadata$Biosample.term.name)


if(file.exists("../data/dataframes_for_all_TF_data.Rda")){
	load("../data/dataframes_for_all_TF_data.Rda")
	}else{
	
#for 1 cell line first
	cell_line<- "CH12.LX"
	TF1<- split$CH12.LX
	TF1<-gsub("-mouse", "", TF1)
	TF1<- unique(TF1)
	#get rid of the ones that you didn't use e.g. no pwm
	data1<-list()
	data1[TF1]<- lapply(TF1, getCellLineData)
	save(data1, file=paste0("../data/dataframe_",cell_line,".Rda"))
	
#second line
	cell_line<- "MEL"
	TF2<-split$MEL
	TF2<-gsub("-mouse", "", TF2)
		TF2<- unique(TF2)

	data2<-list()
	data2[TF2]<- lapply(TF2, getCellLineData)
	save(data2, file=paste0("../data/dataframe_",cell_line,".Rda"))

#3
	cell_line<- "C2C12"
	TF3<- split$C2C12
	TF3<-gsub("-mouse", "", TF3)
	TF3<- unique(TF3)

	data3<-list()
	data3[TF3]<- lapply(TF3, getCellLineData)
	save(data3, file=paste0("../data/dataframe_",cell_line,".Rda"))

#4
	cell_line<- "E14TG2a.4"
	TF4<- split$E14TG2a.4
	TF4<-gsub("-mouse", "", TF4)
	TF4<- unique(TF4)

	data4<-list()
	data4[TF4]<- lapply(TF4, getCellLineData)
	save(data4, file=paste0("../data/dataframe_",cell_line,".Rda"))

#5
	cell_line<- "ES-Bruce4"
	TF5<- split$"ES-Bruce4"
	TF5<-gsub("-mouse", "", TF5)
	TF5<- unique(TF5)

	data5<-list()
	data5[TF5]<- lapply(TF5, getCellLineData)
	save(data5, file=paste0("../data/dataframe_",cell_line,".Rda"))

#6
	cell_line<- "ES-E14"
	TF6<- split$"ES-E14"
	TF6<-gsub("-mouse", "", TF6)
	TF6<- unique(TF6)

	data6<-list()
	data6[TF6]<- lapply(TF6, getCellLineData)
	save(data6, file=paste0("../data/dataframe_",cell_line,".Rda"))

#7
	cell_line<- "G1E"
	TF7<- split$G1E
	TF7<-gsub("-mouse", "", TF7)
	TF7<- unique(TF7)

	data7<-list()
	data7[TF7]<- lapply(TF7, getCellLineData)
	save(data7, file=paste0("../data/dataframe_",cell_line,".Rda"))

#8
	cell_line<- "G1E-ER4"
	TF8<- split$"G1E-ER4"
	TF8<-gsub("-mouse", "", TF8)
	TF8<- unique(TF8)

	data8<-list()
	data8[TF8]<- lapply(TF8, getCellLineData)
	save(data8, file=paste0("../data/dataframe_",cell_line,".Rda"))

data<- list(CH12=c(data1), MEL=c(data2), C2C12=c(data3),E14TG2a.4=c(data4), "ES-Bruce4"=c(data5), "ES-E14"=c(data6), G1E=c(data7), "G1E-ER4"=c(data8) )


save (data, file="../data/set2/dataframes_for_all_TF_data.Rda")
}


#########
#CH12

if(file.exists("../scripts/CH12.LX_optimal_data.csv")){
	read.table("../scripts/CH12.LX_optimal_data.csv", sep=",", header=T)
	}else{
	for (i in 1:length(data1)){
		
		data1[[i]] <- cbind("QDA"=rownames(data1[[i]]), data.frame(data1[[i]], row.names=NULL))

	}

	df1<- data.frame()
	for (i in 1:length(data1)){

		line1<-data1[[i]][which.max(data1[[i]][,2]),]
		df1<- rbind (df1, line1)

	}

rownames(df1)<- TF1
write.csv(df1, "../scripts/CH12.LX_optimal_data.csv", row.names=T)
}
#########
#MEL


if(file.exists("../scripts/MEL_optimal_data.csv")){
	read.table("../scripts/MEL_optimal_data.csv", sep=",", header=T)
}else{
	for (i in 1:length(data2)){
		
		data2[[i]] <- cbind("QDA"=rownames(data2[[i]]), data.frame(data2[[i]], row.names=NULL))

	}

	df2<- data.frame()
	for (i in 1:length(data2)){

		line2<-data2[[i]][which.max(data2[[i]][,2]),]
		df2<- rbind (df2, line2)

	}

rownames(df2)<- TF2
write.csv(df2, paste0("../scripts/MEL_optimal_data.csv"), row.names=T)
}

#C2C12
for (i in 1:length(data3)){
		
		data3[[i]] <- cbind("QDA"=rownames(data3[[i]]), data.frame(data3[[i]], row.names=NULL))

	}

	df3<- data.frame()
	for (i in 1:length(data3)){

		line3<-data3[[i]][which.max(data3[[i]][,2]),]
		df3<- rbind (df3, line3)

	}

rownames(df3)<- TF3
write.csv(df3, paste0("../scripts/C2C12_optimal_data.csv"), row.names=T)

#E14TG2a.4
for (i in 1:length(data4)){
		
		data4[[i]] <- cbind("QDA"=rownames(data4[[i]]), data.frame(data4[[i]], row.names=NULL))

	}

	df4<- data.frame()
	for (i in 1:length(data4)){

		line4<-data4[[i]][which.max(data4[[i]][,2]),]
		df4<- rbind (df4, line4)

	}

rownames(df4)<- TF4
write.csv(df4, paste0("../scripts/E14TG2a.4_optimal_data.csv"), row.names=T)

#ES-Bruce4
for (i in 1:length(data5)){
		
		data5[[i]] <- cbind("QDA"=rownames(data5[[i]]), data.frame(data5[[i]], row.names=NULL))

	}

	df5<- data.frame()
	for (i in 1:length(data5)){

		line5<-data5[[i]][which.max(data5[[i]][,2]),]
		df5<- rbind (df5, line5)

	}

rownames(df5)<- TF5
write.csv(df5, paste0("../scripts/ES-Bruce4_optimal_data.csv"), row.names=T)

#ES-E14
for (i in 1:length(data6)){
		
		data6[[i]] <- cbind("QDA"=rownames(data6[[i]]), data.frame(data6[[i]], row.names=NULL))

	}

	df6<- data.frame()
	for (i in 1:length(data6)){

		line6<-data6[[i]][which.max(data6[[i]][,2]),]
		df6<- rbind (df6, line6)

	}

rownames(df6)<- TF6
write.csv(df6, paste0("../scripts/ES-E14_optimal_data.csv"), row.names=T)

#G1E
for (i in 1:length(data7)){
		
		data7[[i]] <- cbind("QDA"=rownames(data7[[i]]), data.frame(data7[[i]], row.names=NULL))

	}

	df7<- data.frame()
	for (i in 1:length(data7)){

		line7<-data7[[i]][which.max(data7[[i]][,2]),]
		df7<- rbind (df7, line7)

	}

rownames(df7)<- TF7
write.csv(df7, paste0("../scripts/G1E_optimal_data.csv"), row.names=T)

#G1E-ER4
for (i in 1:length(data8)){
		
		data8[[i]] <- cbind("QDA"=rownames(data8[[i]]), data.frame(data8[[i]], row.names=NULL))

	}

	df8<- data.frame()
	for (i in 1:length(data8)){

		line8<-data8[[i]][which.max(data8[[i]][,2]),]
		df8<- rbind (df8, line8)

	}

rownames(df8)<- TF8
write.csv(df8, paste0("../scripts/G1E-ER4_optimal_data.csv"), row.names=T)

#######################
#lambda frequency plot
#listoflambda<- unlist(list(df1[,3],df2[,3],df3[,3],df4[,3],df5[,3],df6[,3],df7[,3],df8[,3]))

data1<-data.frame( cell=c("CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX"), lambda=c(2.50, 1.00, 1.00, 0.75, 2.00, 0.25, 1.25, 1.00, 0.50, 1.00, 
1.50, 0.75, 1.00, 1.75, 2.75, 0.75, 0.75, 2.75, 1.50, 0.75, 1.25, 0.50, 2.25, 0.75, 2.50, 1.75, 1.00, 2.25))

data2<-data.frame( cell=c("MEL", "MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL",
"MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL"), lambda=c(4 ,4,2.5,0.5,2.5,0.5,1.25,0.5,4.75, 0.75,
1.5,1.5, 0.75, 0.75, 0.75, 2.5, 1,1,1,0.5, 1, 1.75, 4,1.5,0.5, 2.25, 0.75, 0.75, 0.5))
data3<- data.frame( cell=c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","C2C12"), lambda=c(1.25, 0.75, 0.75, 1.25, 2.75, 0.75, 1.25))
data4<- data.frame(cell=c("E14TG2a.4","E14TG2a.4","E14TG2a.4"), lambda=c(0.75, 1.75, 1.25))

data5<- data.frame(cell=c("ES-Bruce","ES-Bruce"), lambda=c(1.25, 1.50))
data6<- data.frame(cell=c("ES-E14","ES-E14","ES-E14","ES-E14"), lambda=c( 0.50, 0.50, 0.75, 1.00))
data7<- data.frame(cell=c("G1E","G1E","G1E","G1E"), lambda=c(0.75, 0.5, 0.75, 0.75))
data8<- data.frame(cell=c("G1E-ER4"), lambda=c( 0.75))

lambdadata<-rbind(data1, data2, data3, data4, data5, data6, data7, data8)

library(ggplot2)
ggplot(lambdadata, aes(x=lambda, color=cell))+
geom_histogram(binwidth=0.25,fill="white", alpha=0.5, position="dodge")+
facet_grid(cell ~ .)+
 ggtitle("Frequency of optimal specificity values")+
   theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
	plot.title = element_text(color="black", size=14, face="bold",hjust=0.5))+
ggsave("../imgs/set2/lambdahist_fig3_v4.pdf",width = 10, height = 10)

#######################
#bm frequency plot
data1<-data.frame( cell=c("CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX"), abundance=c( 1e+06, 2e+05, 1e+06, 1e+04, 1e+06, 5e+05, 1e+06,
 5e+03, 1e+05, 1e+05, 1e+06, 2e+05,5e+04 ,5e+05, 1e+06, 2e+04, 5e+04, 1e+06, 1e+05, 1e+04, 2e+05, 1e+05, 1e+06, 5e+04,1e+06, 1e+06, 2e+04, 1e+06))

data2<-data.frame( cell=c("MEL", "MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL",
"MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL"), abundance=c( 1e+06, 1e+06, 1e+06, 10000, 1e+06, 2e+05, 5000,
 10000, 1e+06,50000, 1e+06, 1e+06, 50000, 1e+05, 50000, 1e+06, 1e+05, 20000,20000,5e+05 ,1e+05, 1e+06, 1e+06, 1e+06, 10000, 5e+05, 50000,10000, 2000))

data3<- data.frame( cell=c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","C2C12"),abundance=c(1e+06, 50000, 1e+05, 1e+06, 1e+06, 2e+05, 1e+06))

data4<- data.frame(cell=c("E14TG2a.4","E14TG2a.4","E14TG2a.4"), abundance=c(1e+05, 1e+06, 1e+06))

data5<- data.frame(cell=c("ES-Bruce","ES-Bruce"), abundance=c(1e+05, 1e+06))

data6<- data.frame(cell=c("ES-E14","ES-E14","ES-E14","ES-E14"), abundance=c(1e+04, 1e+06, 1e+06, 5e+04))

data7<- data.frame(cell=c("G1E","G1E","G1E","G1E"), abundance=c(5e+05, 1e+06, 2e+05, 5e+05))

data8<- data.frame(cell=c("G1E-ER4"), abundance=c(1e+06))

abundancedata<-rbind(data1, data2, data3, data4, data5, data6, data7, data8)

ggplot(abundancedata, aes(x=abundance, color=cell))+
geom_histogram(binwidth=0.25,fill="white", alpha=0.5, position="dodge")+
scale_x_log10(breaks=c(1000,2000, 5000, 10000, 20000,50000,1e+04 ,1e+05,2e+05,5e+05,1e+06))+
facet_grid(cell ~ .)+
  ggtitle("Frequency of optimal abundance values")+
   theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
	plot.title = element_text(color="black", size=14, face="bold",hjust=0.5))+
ggsave("../imgs/set2/abundancehist_fig3_v5.pdf",width = 10, height = 10)

#######################
#optimal QDA frequency plot

data1<-data.frame( cell=c("CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX",
"CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX","CH12.LX"), optimalQDA=c( 0.9,  0.95, 0.8,  0.9, 0.7,  0.8,  0.9 , 0.8 , 
0.8,  0.9 , 0.9 , 0.8 , 0.8 , 0.9,  0.9,0.9,  0.7,  0.8 , 0.9,  0.9,  0.8, 0.9 , 0.9 , 0.8,  0.9,  0.95, 0.7 , 0.9))


data2<-data.frame( cell=c("MEL", "MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL",
"MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL","MEL"), optimalQDA=c(0.9 , 0.9, 0.95,0.9 , 0.95 ,0.8,  0.95, 
0.9 , 0.9 , 0.9,  0.9 ,0.9 ,0.95,0.9 ,0.9, 0.95,0.9, 0.7  ,0.9,  0.9 , 0.9 , 0.9, 0.9 , 0.9,  0.7 , 0.95, 0.9,  0.8,  0.9))

data3<- data.frame( cell=c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","C2C12"),optimalQDA=c(0.3, 0 ,  0.3, 0.3, 0.8, 0.7 ,0.7))

data4<- data.frame(cell=c("E14TG2a.4","E14TG2a.4","E14TG2a.4"),optimalQDA=c(0  , 0.3 ,0.1))

data5<- data.frame(cell=c("ES-Bruce","ES-Bruce"),optimalQDA=c(0.7 ,0.8))

data6<- data.frame(cell=c("ES-E14","ES-E14","ES-E14","ES-E14"),optimalQDA=c(0.6, 0.8 ,0.8, 0.9))

data7<- data.frame(cell=c("G1E","G1E","G1E","G1E"),optimalQDA=c( 0, 0.9,  0.99 ,0.9))

data8<- data.frame(cell=c("G1E-ER4"),optimalQDA=c(0.1))
optimalQDAdata<-rbind(data1, data2, data3, data4, data5, data6, data7, data8)

ggplot(optimalQDAdata, aes(x=optimalQDA, color=cell))+
geom_histogram(breaks= c(0, 0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),fill="white", alpha=0.5, position="dodge")+
facet_grid(cell ~ .)+
scale_x_discrete(limit = c(0, 0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99))+
  ggtitle("Frequency of optimal QDA values")+
  theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
	plot.title = element_text(color="black", size=14, face="bold",hjust=0.5))+
ggsave("../imgs/set2/optQDAhist_fig3_v6.pdf",width = 10, height = 10)

