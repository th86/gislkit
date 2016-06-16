#Based on the MATLAB code by Brunet et al., 2004.
#http://www.pnas.org/content/101/12/4164.abstract

nmf<-function( V, num_rank , w=NULL, maxIter=10000, threshold=400, normalize.h=FALSE ) {
	inc=0
	if(length(w)==0)
		w<-matrix( runif( nrow(V)*num_rank,0,10), nrow(V), num_rank)

	h<-matrix( runif( ncol(V)*num_rank,0,1), num_rank ,ncol(V))

	conMat<-matrix(0,ncol(V),ncol(V))
	conMat.prev<-conMat

	if(normalize.h==TRUE){
		h<-h/colSums(h)
	}

	for( i in 1:maxIter ){

		x1<-matrix( t(colSums(w)) ,  num_rank , ncol(V) )
		h<-h*( t(w)%*%(V/(w%*%h)))/x1

		x2<-matrix( t(rowSums(h)) ,  nrow(V) , num_rank )
		w<-w*( (V/(w%*%h))%*%t(h))/x2

		#find the largest in every column of H
		index.max<-apply(h,2,which.max) 

		#conMat<-w%*%h
		#epsilon<-sum( abs( conMat-V) )
		#if( epsilon<threshold)
		#	break

		mat1=matrix( index.max ,ncol(V), ncol(V) )
		mat2=t(matrix( index.max , ncol(V), ncol(V) ))
		conMat<-matrix( as.numeric(mat1==mat2), ncol(V), ncol(V) )

		if( sum(conMat!=conMat.prev)==0) {
			inc=inc+1
		}else{
			inc=0
		}
		if( inc> threshold)
			break;

		conMat.prev<-conMat

	}#end for

	#V.prime<-w%*%h
	
	listOutput<-list(W=w,H=h)
	return( listOutput )
}