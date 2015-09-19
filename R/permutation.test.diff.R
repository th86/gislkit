permutation.test.diff<-function(con_vec,exp_vec,iteration=1000000,verbose=FALSE){
	truth=c(con_vec,exp_vec)
	truth_label=c(rep(0,length(con_vec)), rep(1,length(exp_vec)) )

	names(truth)=truth_label
	ci_truth=equalci( truth , cbind(truth_label,1-truth_label))

	ci_list=rep(0,iteration)
	for(i in 1:iteration){
		ci_list[i]=equalci( sample(truth,length(truth)) , cbind(truth_label,1-truth_label))
		if(i%%1000==1 & verbose==TRUE)
			cat(i,ci_list[i],"\n")
	}

	pval=length(which(ci_list>ci_truth ))/iteration
	if(verbose==TRUE)
		cat("P-value=", pval,"\n")
	
	return( pval )
}