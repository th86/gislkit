#Read a numeric matrix from a csv
read.csv<-function( fn  ){
	ge<-read.table(fn, sep=",", header=T, stringsAsFactors=FALSE)
	rn=make.unique(ge[,1])
	ge<-data.matrix(ge[,-1])
	rownames(ge)=rn
	return(ge)
}

#Read TCGA biotab files
biotab.read<-function( fn  ){
    ge<-data.matrix(read.table(fn,sep="\t", header=T, comment.char=""))
    rn<-as.character(read.table(fn,sep="\t", header=T, comment.char="")[,1])
    #rn=ge[,1]
    ge<-ge[,-1]
    rownames(ge)=rn
    ge<-ge[-1,]
    return(ge)
}
