#Attractor Finding Algorithm for Multiple Seeds
#
# Created by: Taihsien Ouyang (to2232@columbia.edu)
# Reviewed by: Kaiyi Zhu (kz2232@columbia.edu)
# Last Modified 2017 Apr. 19

findMultipleAttractors <- function(dataMatrix, attractorObjects=NULL, a=5, maxIter = 500, epsilon=1E-7, bin = 6, so = 3, NumTopFeature=NumTopFeature, negateMI = TRUE, verbose = TRUE, filterDominatingSeed = TRUE, stopWhenMaxIterReached = FALSE){

	if(is.null(attractorObjects$seedList) == TRUE) stop("Please provide a valid seed list.")

	seedList = attractorObjects$seedList
	attractorSignatureList = attractorObjects$attractorSignatureList
	
	mi <- getAllMIWz(dataMatrix, dataMatrix[seedList[1],], bin=bin, so=so, negateMI = negateMI) #Compute the MI vector
	mi[ mi < 0] = 0
	w <-  mi^a / sum( mi^a )                                                       				#Convert MI to the weight vector
	metagene <- w %*% dataMatrix                                                  				#Compute the metagene using the weight vector
	premi = mi

	#If the score of the second-ranked gene is less 0.5, discard this attractor
  	if(filterDominatingSeed == TRUE & sort(mi,decreasing=TRUE)[2] < 0.5){
  		seedList=seedList[-1]
  		if(verbose==TRUE) print("The attractor is dominated.")
 		return(list(seedList=seedList, attractorSignatureList=attractorSignatureList))
	}

	iter = 1
	while(iter < maxIter){

		mi <- getAllMIWz(dataMatrix, metagene, bin=bin, so=so, negateMI=negateMI)
		mi[ mi < 0] = 0
		w <-  mi^a / sum( mi^a)
		metagene <- w %*% dataMatrix

		delta <- maxAbsDiffTopFeatures(mi, premi, NumTopFeature = NumTopFeature)
		#Print the top 10 features at each iteration
		if(verbose == TRUE){
			cat("Iteration ", (iter+1), "\tDelta = ", delta, "\n", sep="")
			print(mi[order(mi, decreasing=TRUE)[1:10]])
		}

		#Terminate iterations if either of the conditions are valid:
		seedRank <- rank(-mi)[seedList[1]]

		#Call the attractor converged when maximum absolute difference of the top scores between two consecutive iterations < epsilon
		if( is.na(delta) == TRUE ) break

		if( delta < epsilon ){

			#If the attractor is not merged in this iteration, add the attractor to the known attractor list
			newSignatureID = length(attractorSignatureList)+1
			attractorSignatureList[[newSignatureID]] = mi

			if(verbose == TRUE) cat(seedList[1],"is a new attractor\n" )
	
			break

		}else{

			premi <- mi
			iter <- iter + 1

		}

		#Stop iterating if reaches the maximum of iterations.
		if(iter >= maxIter & stopWhenMaxIterReached == TRUE)	stop(paste(seedList[1], "reached the maximum of iterations."))

	}#End of while

	seedList=seedList[-1]

	return (list(seedList=seedList, attractorSignatureList=attractorSignatureList))
}
