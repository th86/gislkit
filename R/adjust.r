#Tools for adjusting the predcitions
#To map the predictions to [0,1] and shift the center by a threhold.
mapScore<-function(pred_scoring, threshold, new_center=0.5){
		pred_scoring=pred_scoring+threshold
		pred_scoring[pred_scoring>0]=pred_scoring[pred_scoring>0]/max(pred_scoring[pred_scoring>0])*0.5 
		pred_scoring[pred_scoring<0]=-1*pred_scoring[pred_scoring<0]/min(pred_scoring[pred_scoring<0])*0.5
		pred_scoring=pred_scoring+new_center
	return(pred_scoring)
}

#To find threhsold for the best BAC
findThreshold<-function( prediction, response, sweep_max=2, resolution=10000, doPlot=FALSE   ){
	pred_interval=(1:(sweep_max*resolution))/resolution
		bac_roc<-rep(0,  sweep_max*resolution )
		for(i in 1:length(pred_interval)   ){
			bac_roc[i]= bac(prediction, response ,threshold= pred_interval[i] )
		}
	
	if(doPlot==TRUE){
		plot(pred_interval , bac_roc, type="l", col=4, main="Threshold vs. BAC", xlab="Threshold",ylab="BAC" )
		points(pred_interval[ which.max(bac_roc) ] ,max(bac_roc), cex=0.5,col=2,pch=20)
	}

	thresholdSet=c(which(bac_roc>(max(bac_roc)-0.001)))
	return(list( score=max(bac_roc), thresholdSet=thresholdSet  ) )
}

#To fit the predicted ranks to the known survival information
fitTime<-function( prediction, prediction_known ,known_survival,halfwindow=5, doPlot=FALSE, increment=0.1  ){
 idx = which(known_survival[, 2] == 1)
    prediction_d = prediction_known[idx]
    known_survival_d = known_survival[idx, ]
    names(prediction_d) = known_survival_d[, 1]
    prediction_d_sorted = sort(prediction_d)
    smooth_table <- rep(0, length(prediction))
    for (i in 1:length(prediction)) {
        larger = as.numeric(names(which(prediction_d_sorted < 
            prediction[i])))
        if (length(larger) > halfwindow) 
            larger = larger[(length(larger) - (halfwindow - 1)):length(larger)]
        smaller = as.numeric(names(which(prediction_d_sorted > 
            prediction[i])))
        if (length(smaller) < halfwindow) {
            smaller = smaller[1:length(smaller)]
        }
        else {
            smaller = smaller[1:halfwindow]
        }
        cat(i, "length", length(larger), length(smaller), "\n")
        smooth_table[i] = mean(c(smaller, larger), na.rm = TRUE)
    }
    pred_scoring_rank = rank(prediction)
    pred_scoring_sorted = sort(prediction)
    smooth_table_sorted = sort(smooth_table)



    predicted_survival = smooth_table_sorted[pred_scoring_rank]


    for (i in 1:length(predicted_survival)) 
    	if (i > 1) 
    		if( length(which(predicted_survival==predicted_survival[i]))>1)
    			predicted_survival[which(predicted_survival==predicted_survival[i])]= predicted_survival[which(predicted_survival==predicted_survival[i])]+increment*(1:length(which(predicted_survival==predicted_survival[i])))

    if (doPlot == TRUE) {
        plot(prediction, predicted_survival, main = "Fitted Survival", 
            cex.main = 0.7, xlab = "Predicted Score", ylab = "Survival", 
            col = 4, pch = 20, cex = 0.6)
    }
    return(predicted_survival)
}

