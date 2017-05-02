loglik <-
function(object){
  r<-object$r
  z<-object$ord.z
  p<-object$p
  Rt<-object$Rt

  if(p==0) Rt0<-Rt

  if(p>0){
    if(p>1) Rt0<-Rt*exp(-sum(object$coefficients*apply(z,2,mean)))
    if(p==1) Rt0<-Rt*exp(-sum(object$coefficients*mean(z)))
  }

  if(p==0) ezb<-1
  if(p>0) ezb<-as.numeric(exp(as.matrix(z)%*%as.matrix(object$coefficients)))

  rt0<-diff(c(0,Rt0))
  if(r>0) rt0<-rt0/(1+r*Rt0*ezb)

  if(r==0) lls<--sum(Rt0*ezb)
  if(r>0) lls<-sum(-log(1+r*Rt0*ezb)/r)
  loglik<-sum(object$ord.delta*log(rt0*ezb*object$ord.delta+(1-object$ord.delta)))+lls
return(loglik)
}
