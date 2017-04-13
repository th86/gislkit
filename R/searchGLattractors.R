#Genomically Localized Attractor Search Algorithm
#
# Created by: Taihsien Ouyang (to2232@columbia.edu)
# Reviewed by: Kaiyi Zhu (kz2232@columbia.edu)
# Last Modified 2017 Apr. 13

#Arguments
#ge.raw:	Input matrix. The rownames must be gene symbols
#maxIter:	Maximum of iterations allowed for the attractor finding process
#epsilon:	Threshold for determining if an attractor is converged
#NumTopFeature:	Number of the features for computing the maximum absolute difference between two vectors
#seedRankThreshold:	Threshold of the seed's rank for determining if an attractor should be terminated
#strengthRange:	The rank of the sorted score that is used as the strength of an attractor
#windowSize:	Number of top features that are printed in the summary report
#alpha:	exponent for genomically localized attractors

#data("grch37.geneymbol") 
#grch37_genesymbol

searchGLattractors<-function(ge.raw, genome, maxIter=500, epsilon=1e-7, NumTopFeature=20, outputLength=50, windowSize=50, alpha=2, seedRankThreshold=1, strengthRange=2){

	print(ls())
	#Prepare the reference genome
	colnames(genome)[2] = "Chr"	
	genomeLocation = 1:nrow(genome)
	names(genomeLocation) = rownames(genome)

	#Sort the genes by their genomic location
	rownames_sorted = intersect ( rownames( genome ) ,rownames(ge.raw) )
	ge = ge.raw[rownames_sorted,] 
	rm(ge.raw)
	geneLocal = rownames(ge)

	#Exhaustive search for the GL attractors in ge, using a=2
	cat("First-pass search\n")
	attractorList=list(seedList=geneLocal, attractorSignatureList=list())		
	while(length(attractorList$seedList) > 0){
		attractorList <- findMultipleGLAttractors(ge, attractorList, refGenome=genome, a=alpha, windowSize=windowSize, maxIter=maxIter, epsilon=epsilon, bin=6, so=3, NumTopFeature=NumTopFeature, negateMI=TRUE, verbose=FALSE)
		cat("Remaining seeds:",length(attractorList$seedList),"\tNumber of attractors:" , length(attractorList$attractorSignatureList) ,"\n" )
	}

	#Convert the valid attractors to pre-neighborhoods
	attractorTopList = 1:length(attractorList$attractorSignatureList) #The indices of the attractors
	preneighborhoodGeneList = list()
	preneighborhoodFlaggedGeneList = list()
	preneighborhoodNumber = 0

	while(length(attractorTopList) > 0){ #does not work if there are no attractors

		strength = sort(attractorList$attractorSignatureList[[attractorTopList[1]]], decreasing=TRUE)[strengthRange]

		if(strength > 0.5){ #2017 Jan. 10
			preneighborhoodNumber = preneighborhoodNumber+1
			preneighborhoodStart <- max(1,genomeLocation[names(attractorList$attractorSignatureList)[attractorTopList[1]]]-floor(windowSize/2))
			preneighborhoodEnd <- min(genomeLocation[names(attractorList$attractorSignatureList)[attractorTopList[1]]]+floor(windowSize/2), nrow(genome) )
			preneighborhoodGeneList[[preneighborhoodNumber]] = rownames(genome)[preneighborhoodStart:preneighborhoodEnd]
			preneighborhoodFlaggedGeneList[[preneighborhoodNumber]] <- names( which(attractorList$attractorSignatureList[[attractorTopList[1]]] > 0.5)  )
		}

		attractorTopList <- setdiff(attractorTopList, attractorTopList[1])

	}#end of while

	if(is.null(preneighborhoodFlaggedGeneList) == TRUE)
		stop("No valid pre-neighborhoods.")	
	

	#Detecting the preneighborhood with co-expression
	preneighborhoodTopList = 1:length(preneighborhoodGeneList) #The indices of pre-neighborhoods
	GeneList = list()			#Neighborhoods
	FlaggedGeneList = list() 	#Flagged features of the neighborhoods

	#Use the first valid pre-neighborhood as a neighborhood
	numberNeighborhoods = 1
	GeneList[[numberNeighborhoods]] = preneighborhoodGeneList[[1]]
	FlaggedGeneList[[numberNeighborhoods]] = preneighborhoodFlaggedGeneList[[1]]
	preneighborhoodTopList <- setdiff(preneighborhoodTopList, preneighborhoodTopList[1])

	GeneList.this <- GeneList[[numberNeighborhoods]]
	FlaggedGeneList.this <- FlaggedGeneList[[numberNeighborhoods]] 

	#Merge the pre-neighborhoods (attractorList) to the neighborhoods defined by GeneList and FlaggedGeneList
	cat("Merge the pre-neighborhoods\n")
	while( length(preneighborhoodTopList) > 0  ){

		dropList=NULL
		for(itr in preneighborhoodTopList){

			if( length( intersect(FlaggedGeneList.this, preneighborhoodFlaggedGeneList[[itr]] ) ) > 0 ){
				GeneList.this <- union(GeneList.this, preneighborhoodGeneList[[itr]])
				FlaggedGeneList.this <- union(FlaggedGeneList.this, preneighborhoodFlaggedGeneList[[itr]])
				dropList <- union(dropList,itr) #Tag the pre-neighborhoods that are overlapped with neighborhood.this , which will be dropped later
				cat("pre-neighborhood", itr , "is merged to", "neighborhood", numberNeighborhoods, ".\n")
			}

		}#end of for

		#Remove the scanned pre-neighborhoods that are overlapped with neighborhood.this 
		preneighborhoodTopList <- setdiff( preneighborhoodTopList, dropList ) 

		#Use the next un-merged pre-neighborhood as the next neighborhood.this, when no pre-neighborhoods can be merged to any variants of neighborhood.this.
		if(	is.null(dropList) == TRUE ){

			#Update this neighborhood
			GeneList[[numberNeighborhoods]] = GeneList.this
			FlaggedGeneList[[numberNeighborhoods]] = FlaggedGeneList.this
			#Initialize the next neighborhood
			numberNeighborhoods = numberNeighborhoods+1
			GeneList[[numberNeighborhoods]] = preneighborhoodGeneList[[preneighborhoodTopList[1]]]
			FlaggedGeneList[[numberNeighborhoods]] = preneighborhoodFlaggedGeneList[[preneighborhoodTopList[1]]]
			#Reset the register of the neighborhood 
			GeneList.this<-GeneList[[numberNeighborhoods]]
			FlaggedGeneList.this = FlaggedGeneList[[numberNeighborhoods]] 
			#Update the list of remaining pre-neighborhoods
			preneighborhoodTopList <- setdiff( preneighborhoodTopList, preneighborhoodTopList[1] ) 

		}

	}#end of while

	#Exhaustive search using the neighborhoods of the co-expressed genes
	cat("Second-pass search\n")
	attractorOutputList=list()
	for(GeneList_iter in 1:length(GeneList)){

		if( is.null(GeneList[[GeneList_iter]]) == FALSE  ){

			seedListLocal = intersect(GeneList[[GeneList_iter]],rownames(ge))
			attractorListExt = list(seedList=seedListLocal, attractorSignatureList=list())
			cat("Neighborhood",GeneList_iter,":", attractorListExt$seedList[1:min(5, attractorListExt$seedList )] ,"\n")
		  	
		  	while(length(attractorListExt$seedList) > 1 ){
			   attractorListExt <- findMultipleAttractors(ge[seedListLocal,], attractorListExt,  a=alpha, maxIter=maxIter, epsilon=epsilon, bin=6, so=3, NumTopFeature=NumTopFeature, negateMI=TRUE, seedRankThreshold = length(seedListLocal), verbose=FALSE)
			    cat("Remaining seeds:",length(attractorListExt$seedList),"\tNumber of attractors:" , length(attractorListExt$attractorSignatureList) ,"\n" )
		    }

		    if(length(attractorListExt$attractorSignatureList)>0){
		      attractorStrength = rep(0,length(attractorListExt$attractorSignatureList))
		      for(topgene_iter in 1:length(attractorListExt$attractorSignatureList) )
		      	attractorStrength[topgene_iter] <- sort(attractorListExt$attractorSignatureList[[topgene_iter]],decreasing=TRUE)[5]

		      cat("Found", names(sort(attractorListExt$attractorSignatureList[[which.max(attractorStrength)]],decreasing=TRUE))[1:5],"\n"  )
		      attractorOutputList[[GeneList_iter]]=attractorListExt$attractorSignatureList[[which.max(attractorStrength)]]
		    }

		}

	} #end of for

	resultObject<-list(attractorList=attractorOutputList, GeneList=GeneList, FlaggedGeneList=FlaggedGeneList)

	return(resultObject)
}
