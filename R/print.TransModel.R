print.TransModel <-
function(x,...){
	cat("Call:\n")
	print(x$call)
  if(x$p==0){ 
    cat("No covariates/ null model.","\n")
    cat("Loglik: ",round(loglik(x),3),".\n")
  }
  if(x$p>0){
	cat("\nCoefficients:\n")
	print(x$coefficients)
	cat("\nCovariance Matrix:\n")
	print(x$vcov)
      cat("Loglik: ",round(loglik(x),3),", AIC=",round(-2*loglik(x)+2*x$p,3),".\n")
  }
	class(x)<-"TransModel"
}
