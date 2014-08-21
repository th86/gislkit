direction.correction<-function(pred_in, ref_in, rev=-1){
	
	if(cor(pred_in, ref_in)*rev<0){
		pred_out=-1*pred_in
		cat("Reversed prediction. Corrected.\n")
	}else{
		pred_out=pred_in
		cat("Correct prediction.\n")
	}

	return(pred_out)
}