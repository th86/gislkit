#Compute the maximum absolute difference of the top genes of two vectors
maxAbsDiffTopFeatures<-function(vectorX, vectorY, NumTopFeature=20){

	topFeatureList=intersect(names(sort(vectorX, decreasing=TRUE)), names(sort(vectorY, decreasing=TRUE)) )[1:NumTopFeature]

	return( max(abs(vectorX[topFeatureList] - vectorY[topFeatureList])) )
}

#Compute the maximum absolute difference of the top genes of a vector and all the vectors in a list
allMaxAbsDiffTopFeatures <- function(vectorX, signatureList, NumTopFeature=20){
	maxabsScoreDiffVector <- rep(1,length(signatureList))

	for(i in 1:length(signatureList)){
		maxabsScoreDiffVector[i] <- maxAbsDiffTopFeatures(vectorX, signatureList[[i]], NumTopFeature)
	}

	return(maxabsScoreDiffVector)
}
