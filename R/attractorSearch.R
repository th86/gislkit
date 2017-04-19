attractorSearch<-function(dataMatrix, geneList , a = 5, maxIter = 500, epsilon = 1e-7, bin = 6, so = 3, NumTopFeature=20, negateMI = TRUE, verbose = FALSE){

	seedList<-intersect( rownames(dataMatrix), geneList )

	if( length(seedList) > 0 ){

		cat("Seed list:", seedList, "\n", sep = "\t" )

		attractorList = list(seedList = seedList, attractorSignatureList = list())
		start.time <- Sys.time()
		while(length(attractorList$seedList) > 0){

			attractorList <- findMultipleAttractors(dataMatrix, attractorList, a = a, 
			maxIter = maxIter, epsilon = epsilon, bin = bin, so = so, NumTopFeature = NumTopFeature, negateMI = negateMI, 
			filterDominatingSeed = FALSE, verbose = verbose)
		
		}
		end.time <- Sys.time()
		time.taken <- end.time - start.time
		print(time.taken)
		return(attractorList$attractorSignatureList)

	}else{

		stop("Please provide a valid list of gene name(s)")

	}
}

