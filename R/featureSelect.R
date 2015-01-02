#Randomized Feature Selection
FeatureSelect<-function(covar,time,status,featureSet_next, maxIter=100000){

	featureSet_next=intersect(featureSet_next, colnames(covar))

	score_next=0
	score_prev=0
	featureSet_prev=featureSet_next

	df<-data.frame(cbind(time,status, covar[,featureSet_next] ) )
	colnames(df)[1:2]=c("time","cens")

	score_prev=cv.cox(df,useFULL=TRUE)

	for(i in 1:maxIter){
		featureSet_retained<-sample(featureSet_next, length(featureSet_next)-1)
		featureSet_sampled<-sample(setdiff(colnames(aml),featureSet_next),1 )
		featureSet_next=c(featureSet_retained, featureSet_sampled )

		df<-data.frame(cbind(time, status, covar[,featureSet_next] ) )
		colnames(df)[1:2]=c("time","cens")

			score_next=cv.cox(df,useFULL=TRUE, iteration=1)
				if(is.nan( score_next )==FALSE )
					if(score_next>score_prev){
						featureSet_prev=featureSet_next
						score_prev=score_next
						cat( score_prev, score_prev, "\n" )
						cat( featureSet_next, "\n" )
			}else{
				featureSet_next=featureSet_prev
			}

		if(i%%100==1)
		cat("Iteration:",i,"Feature:",sort(featureSet_prev),"\n")

	}

	cat("Selection:", score_prev,sort(featureSet_prev), "\n" )

    return(c(score_prev,sort(featureSet_prev)))
}

#Batch feature selection
batchFeatureSelect<-function(covar,maxfeatureNumber=10,fileName="featureSelectionList.rda" ){
	featureList=list()
	for(i in 2:maxfeatureNumber)
    	featureList[[i]]=FeatureSelect(covar, c(sample(colnames(aml.imputed),i)))

	save(featureList, file=fileName)
}

#Summarize batch feature selection
summarizeFeatureSelect<-function(pattern=".rda",maxfeatureNumber=10){

	fnList=dir(pattern=pattern)

	maxfeatureList=list()

	for(i in 2:maxfeatureNumber)
		maxfeatureList[[i]]=c(0,0)

	for(fn in fnList){
		load(fn)
		cat(fn,featureList[[maxfeatureNumber]][1],"\n")
		for(i in 2:maxfeatureNumber){
			if(featureList[[i]][1]>maxfeatureList[[i]][1] )
				maxfeatureList[[i]]=c(featureList[[i]])
		}
	}

	return(maxfeatureList)
}


#Coefficient Optimization
coefOptimize<-function(df,time,status,sd=1,maxIter=1000000,coef_prev=NULL){

	if(is.null(coef_prev)==TRUE)
			coef_prev=rnorm(ncol(df))
	
	ci_prev=equalci(aml[,featureSet_next]%*%coef_prev, Surv(time,status))
	coef_next=coef_prev
	for(i in 1:maxIter){

	    coef_select=sample(ncol(df),1)

	    #if(sample(2,1)>1){
	    #    coef_next[coef_select]=coef_prev[coef_select]+0.000001
	    #}else{
	    #    coef_next[coef_select]=coef_prev[coef_select]-0.000001
	    #}

	    coef_next[coef_select]=coef_prev[coef_select] +rnorm(1,mean=0,sd=sd)

	    coef_next=coef_next/max(abs(coef_next))

	    pred=aml[,featureSet_next]%*%coef_next
	    ci_next=equalci(pred, Surv(time,status))
	    cat(i, ci_prev, ci_next, coef_select, coef_next,"\n")

	    if(ci_next>=ci_prev){
	        coef_prev=coef_next
	        ci_prev=ci_next
	        cat("UPDATE",ci_next,"\n")
	    }

	}

	names(coef_prev)=featureSet_next
	cat(ci_prev,coef_prev ,"\n")

	return(coef_prev)
}
