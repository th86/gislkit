#Read a numeric matrix from a csv
readcsvmatrix<-function( fn  ){
	ge<-read.table(fn,sep=",", header=T)
	rn=ge[,1]
	ge<-data.matrix(ge[,-1])
	rownames(ge)=rn
	return(ge)
}