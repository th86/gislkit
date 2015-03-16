#Print top features in each column
topfeat<-function(obj,list_length=10){
	w<-basis(obj)
	for(i in 1:ncol(w) ) cat(i,"\n",names(sort(w[,i],decreasing=T))[1:list_length],"\n",round(sort(w[,i],decreasing=T)[1:list_length],4),"\n")
}
