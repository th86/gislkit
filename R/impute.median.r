#Find the balanced accuracy for a vector of predictions
impute.median=function(mat){
   mat<-as.numeric(as.character(mat)) 
   mat[which(is.na(mat==TRUE))] =median(mat, na.rm=TRUE) 
   return(mat) 
}