print.summary.TransModel <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	if(x$p==0){
    	  cat("No covariates/ null model: ","\n")
    	  cat("n = ",length(x$ord.delta),"\n")
    	  cat("n.event = ",sum(x$ord.delta),"\n")
	  cat("Loglik: ",round(loglik(x),3),".\n")
	}
	if(x$p>0){ 
	  printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
	  cat("Loglik: ",x$loglik,",AIC=",x$AIC,".\n")
	}
	class(x)<-"summary.TransModel"
}
