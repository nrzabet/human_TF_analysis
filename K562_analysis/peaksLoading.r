peaksLoading<-function(x){
    # Laoding files based on files extension

    if(grepl(x=x,pattern=".bed")& !grepl(x=x,pattern=".gff")){

          x<-read.table(x, stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }

    }else if(grepl(x=x, pattern=".Rda")){
        x<-get(load(x))
    }else if(grepl(x=x, pattern=".gff3")){
          x<-read.table(x, skip=30,stringsAsFactors=F)
          if(length(grep(x=x[,1], pattern="chr"))==0){
              x[,1]<-paste0("chr",x[,1])
              x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
          }else{
              x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
          }

    } else if(grepl(x=x, pattern=".gff")){
        x<-read.table(x, skip=30,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }else{
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,4],x[,5]))
        }
    } else{
        x<-read.table(x, header=T,stringsAsFactors=F)
        if(length(grep(x=x[,1], pattern="chr"))==0){
            x[,1]<-paste0("chr",x[,1])
            x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }else{
             x<-GRanges(seqnames=as.character(x[,1]),range=IRanges(x[,2],x[,3]))
        }
    }
    return(x)
}
