#Find the balanced accuracy for a vector of predictions
bac<-function(prediction,observation,threshold=0.5){
	P_set<-length(which(observation==1))
	N_set<-length(which(observation==0))
	TP_set<-length(which( prediction>=threshold & observation==1  ))
	TN_set<-length(which( prediction<threshold & observation==0  ))
	bac<-0.5*TP_set/P_set + 0.5*TN_set/N_set
	return(bac)
}
