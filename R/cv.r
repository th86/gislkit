#cross validated feature selection of the machine learning algorithms
require("e1071")
require("survival")
library("gbm")
#require("glmnet")

#cv.logistic
#cv.logistic<-cv.glmnet(X, resp, nfolds=10, family="binomial" ,type.measure="auc")
#logistic_model<-glm(resp~.,family=binomial("logit"),data=data.frame(X))

cv.svm <- function(data, iterations=100, useFULL=FALSE, fold=0.75, kernel = "linear",  cost=1) {
	ci_list=rep(0,iterations)
	for(iter in 1:iterations){

		if(useFULL==FALSE){
	    	trainSet=data[sample(rownames(data), floor(nrow(data)*fold)),]
	    	testSet=data[setdiff(rownames(data),rownames(trainSet)),]
	 	}else{
	 		trainSet=data
	 		testSet=data	
	 	}

	svm_model<-svm(   resp~., data=trainSet , kernel=kernel, type="C",scale=TRUE , cost=cost ,cross=10, probability=TRUE) #gamma 0.058
    pred<-predict(svm_model, testSet)
    ci_list[iter]<-exactci(pred, Surv(1-testSet[,"resp"], testSet[,"resp"]))

    }
    return(mean(ci_list,na.rm=TRUE) )
}


cv.cox <- function(data, iterations=100, useFULL=FALSE, fold=0.75) {
	ci_list=rep(0,iterations)
	for(iter in 1:iterations){

		if(useFULL==FALSE){
	    	trainSet=data[sample(rownames(data), floor(nrow(data)*fold)),]
	    	testSet=data[setdiff(rownames(data),rownames(trainSet)),]
	 	}else{
	 		trainSet=data
	 		testSet=data	
	 	}

	    cox_model<-tryCatch({
	     				coxph(Surv(time, cens) ~., data = trainSet, control = coxph.control(iter.max = 100))
	     				},warning=function(w){
	     				},error=function(e){
	     					print(e)
	     				},finally= {
	     				})
	    if(is.null(cox_model)==FALSE){
		    cox_pred<-predict(cox_model, testSet)
	    	pred=cox_pred
	     	ci_list[iter]<-exactci(pred, Surv( testSet[,"time"],testSet[,"cens"]) )
	    }
	}
    return(mean(ci_list,na.rm=TRUE) )
}


cv.cox.equal<-function (data, iterations = 100, useFULL = FALSE, fold = 0.75) 
{
    ci_list = rep(0, iterations)
    for (iter in 1:iterations) {
        if (useFULL == FALSE) {
            trainSet = data[sample(rownames(data), floor(nrow(data) * 
                fold)), ]
            testSet = data[setdiff(rownames(data), rownames(trainSet)), 
                ]
        }
        else {
            trainSet = data
            testSet = data
        }
        cox_model <- tryCatch({
            coxph(Surv(time, cens) ~ ., data = trainSet, control = coxph.control(iter.max = 100))
        }, warning = function(w) {
        }, error = function(e) {
            print(e)
        }, finally = {
        })
        if (is.null(cox_model) == FALSE) {
            cox_pred <- predict(cox_model, testSet)
            pred = cox_pred
            ci_list[iter] <- equalci(pred, Surv(testSet[, "time"], 
                testSet[, "cens"]))
        }
    }
    return(mean(ci_list, na.rm = TRUE))
}


cv.gbm <- function(data, iterations=10, useFULL=FALSE, fold=0.75, cv.folds=3,shrinkage=0.05) {
	ci_list=rep(0,iterations)
	for(iter in 1:iterations){

		if(useFULL==TRUE){
	    	trainSet=data[sample(rownames(data), floor(nrow(data)*fold)),]
	    	testSet=data[setdiff(rownames(data),rownames(trainSet)),]
	 	}else{
	 		trainSet=data
	 		testSet=data	
	 	}

	 gbm_model <- gbm(Surv(time, cens) ~., data = trainSet, distribution="coxph" ,cv.folds=cv.folds, shrinkage=shrinkage, verbose = FALSE)
  
     gbm_pred<-predict(gbm_model, testSet)
     pred=gbm_pred
     ci_list[iter]<-exactci(pred, Surv( testSet[,"time"],testSet[,"cens"]) )
    }
    return(mean(ci_list,na.rm=TRUE) )
}
