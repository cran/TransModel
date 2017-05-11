print.TransModel <-
function(x,...){
	cat("Call:\n")
	print(x$call)
  if(x$p==0){ 
    cat("No covariates/ null model.","\n")
  }
  if(x$p>0){
	cat("\nCoefficients:\n")
	print(x$coefficients)
	cat("\nCovariance Matrix:\n")
	print(x$vcov)
  }
	class(x)<-"TransModel"
}
