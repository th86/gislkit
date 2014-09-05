#Randomized Feature Selection
fs<-function( e , remission, featureSet_next, maxIter=10000   ){

	score_full_next=0
	score_full_prev=0
	score_cv_next=0
	score_cv_prev=0

	featureSet_prev=featureSet_next


	    df<-data.frame(cbind(remission[,1],remission[,2] , e[,featureSet_next] ) )
	    colnames(df)[1:2]=c("time","cens")
	    score_full_next=cv.cox.equal(df,useFULL=TRUE)

	for(i in 1:maxIter){

		if(is.nan( score_full_next )==FALSE )
			if(score_full_next>score_full_prev){
				
						featureSet_prev=featureSet_next
						score_full_prev=score_full_next
						score_cv_prev=score_cv_next
						cat(  "S2 ",score_full_prev, score_cv_prev, "\n" ); flush.console()
						cat( featureSet_next, "\n" )
			}else{
				featureSet_next=featureSet_prev
			}

	    featureSet_retained<-sample(featureSet_next, length(featureSet_next)-1)
	    featureSet_sampled<-sample(setdiff(colnames(e),featureSet_next),1 )
	    featureSet_next=c(featureSet_retained, featureSet_sampled )

	    df<-data.frame(cbind(remission[,1],remission[,2] , e[,featureSet_next] ) )
	    colnames(df)[1:2]=c("time","cens")
	    score_full_next=cv.cox.equal(df,useFULL=TRUE)

		if(i%%100==1)
		print(i)

	}

	cat( score_full_prev, score_cv_prev, "\n" )
	cat( featureSet_next, "\n" )

		selectedFeature=list(score=score_full_prev, feature=featureSet_prev)
		
return( selectedFeature )
}