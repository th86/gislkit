#Read a numeric matrix from a csv
readcsvmatrix<-function( fn  ){
	ge<-data.matrix(read.table(fn,sep=",", header=T))
	rn=ge[,1]
	ge<-ge[,-1]
	rownames(ge)=rn
	return(ge)
}