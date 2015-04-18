#Shuffle a matrix by label or by each feature
matrix.shuffle<-function(covar,bylabel=FALSE){

	M<-matrix(0,nrow(covar), ncol(covar))
	rownames(M)<-rownames(covar)
	colnames(M)<-colnames(covar)
	if(bylabel==TRUE){
		M<-M[sample.int(nrow(M)),]
		
	}else{
		for( protein.itr in 1:ncol(covar))
			M[,protein.itr]<-covar[sample.int(nrow(M)), protein.itr]	    
	}

	return(M)
}