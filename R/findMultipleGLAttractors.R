#Genomic Localized Attractor Finding Algorithm With Variable Exponent Using Multiple Seeds
#
# Created by: Taihsien Ouyang (to2232@columbia.edu)
# Reviewed by: Kaiyi Zhu (kz2232@columbia.edu)
# Last Modified 2017 Apr. 13

findMultipleGLAttractors <- function(dataMatrix=NULL, attractorObjects=NULL, refGenome=NULL, windowSize=50, a=6, strengthRank=5, maxIter = 500, epsilon=1E-7, bin = 6, so = 3, NumTopFeature = 20, negateMI = TRUE, verbose = TRUE, filterDominatingSeed = FALSE, stopWhenMaxIterReached = FALSE){

	if(is.null(attractorObjects$seedList) == TRUE) stop("Please provide a valid seed list.")
	if(is.null(refGenome) == TRUE) stop("Please provide a valid reference genome list.")

	seedList <- attractorObjects$seedList
	attractorSignatureList <- attractorObjects$attractorSignatureList

	#Find the genes in the neighborhood of the seed. the range is 2*windowSpan+1
	windowSpan <- floor(windowSize/2)
	centerGeneLoc <- which( rownames(refGenome) == seedList[1] )
	centerGeneChr <- refGenome[centerGeneLoc, "Chr"]
	refGenomeChrList <- rownames(refGenome)[which(refGenome[,"Chr"] == centerGeneChr)]			#Find the window in the chromosome which the seed belongs to
	centerGeneChrLoc <- which( refGenomeChrList == seedList[1] )

	if(length(centerGeneChrLoc) == 0){
		seedList = seedList[-1]
  		if(verbose == TRUE) print("The gene cannot be located.")
 		return(list(seedList=seedList, attractorSignatureList=attractorSignatureList))

	}else{

		if(verbose == TRUE) cat(seedList[1], "is gene #", centerGeneChrLoc, "in", centerGeneChr,"\n")
	}

	geneWindowList <- refGenomeChrList[max(1, centerGeneChrLoc - windowSpan) : min(length(refGenomeChrList), centerGeneChrLoc + windowSpan)]
	dataMatrix <- dataMatrix[intersect(geneWindowList,rownames(dataMatrix)),] 					#Extract the rows by the genomic location

	attractorStrengthMax=0
	attractorWeight=NULL

	if(seedList[1]%in%rownames(dataMatrix) == TRUE){

		a_stepwise = a
			#First Iteration of attractor finding algorithm
			mi <- getAllMIWz(dataMatrix, dataMatrix[seedList[1],], bin=bin, so=so, negateMI = negateMI) #Compute the MI vector
			mi[ mi < 0] = 0
			if(sum(mi) > 0){

				w <-  mi^a_stepwise / sum( mi^a_stepwise )                                                       				#Convert MI to the weight vector
				metagene <- w %*% dataMatrix                                                  				#Compute the metagene using the weight vector
				premi = mi

			}else{ #Jan. 10, 2017 Handling the all zeros

				seedList=seedList[-1]
				if(verbose==TRUE) print("The attractor is dominated.")
		 			return(list(seedList=seedList, attractorSignatureList=attractorSignatureList))

			}

			##If the score of the second-ranked gene is less 0.5, discard this attractor
		  	if(filterDominatingSeed == TRUE & sort(mi,decreasing=TRUE)[2] < 0.5){

		  		seedList = seedList[-1]
		  		if(verbose == TRUE) print("The attractor is dominated.")
		 		return(list(seedList=seedList, attractorSignatureList=attractorSignatureList))
			}

			iter = 1
			while(iter < maxIter){

				mi <- getAllMIWz(dataMatrix, metagene, bin=bin, so=so, negateMI = negateMI)
				mi[ mi < 0] = 0
				w <-  mi^a_stepwise / sum( mi^a_stepwise)
				metagene <- w %*% dataMatrix

				delta <- maxAbsDiffTopFeatures(mi, premi, NumTopFeature=NumTopFeature)
				if(is.na(delta) == TRUE){
		  			if(verbose == TRUE) print("Invalid MI")
		 			break
				}

				#Print the top 10 features at each iteration
				if(verbose == TRUE){
					cat("Exponent ", a_stepwise,"\tIteration ", (iter+1), "\tDelta = ", delta, "\n", sep="")
					print(mi[order(mi, decreasing=TRUE)[1:10]])
				}

				#Terminate iterations if either of the conditions are valid:
				seedRank <-rank(-mi)[seedList[1]]

				#Call the attractor converged when maximum absolute difference of the top scores between two consecutive iterations < epsilon
				if( delta < epsilon ){
					#Remove the attractors which converges within 10 iterations but the seed's rank is larger than the threshold

						attractorWeight <- mi
						if(verbose == TRUE)	cat(seedList[1], "is a new attractor with a =", a_stepwise, "\n" )

					break

				}else{

					premi = mi
					iter = iter + 1

				}

				#Stop iterating if reaches the maximum of iterations.
				if(iter >= maxIter & stopWhenMaxIterReached == TRUE)	stop(paste(seedList[1], "reached the maximum of iterations."))

			}#End of while

		}

	if(is.null(attractorWeight) == FALSE){
		attractorSignatureList[[length(attractorSignatureList)+1]] = attractorWeight
		names(attractorSignatureList)[length(attractorSignatureList)] = seedList[1] #Keep track of the seeds for reconstructing the neighborhood
	}

	seedList = seedList[-1]

	return (list(seedList=seedList, attractorSignatureList=attractorSignatureList))
}
