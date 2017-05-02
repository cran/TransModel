predict.TransModel <-
function(object,...){
  arg<-list(...)
  r<-object$r
  z<-object$z
  if(is.null(arg$newdata)) arg$newdata<-rep(0,object$p)
  arg$newdata<-as.numeric(arg$newdata)
  Rt<-object$Rt
  if(is.null(arg$alpha)) arg$alpha<-0.05

  if(is.null(arg$new.time)) arg$new.time<-object$ord.time[object$ord.delta==1]
  n<-length(arg$new.time)
  p<-object$p

  if(p==0) err<-Rt

  if(p>0){
    if(p>1) Rt0<-Rt*exp(-sum(object$coefficients*apply(z,2,mean)))
    if(p==1) Rt0<-Rt*exp(-sum(object$coefficients*mean(z)))
    err<-Rt0*exp(sum(arg$newdata*object$coefficients))
  }

  if(r==0) St<-exp(-err)
  if(r>0) St<-exp(-1/r*log(1+r*err))
  St.fun<-stepfun(object$ord.time,c(1,St))
  new.St<-St.fun(arg$new.time)
  pred<-data.frame(time=arg$new.time,survival=new.St)

  if(r==0) new.err<-log(-log(new.St))
  if(r>0) new.err<-log((exp(-r*log(new.St))-1)/r)

  if(object$CICB.st){
    new.Ht<-matrix(ncol=n,nrow=nrow(object$rep.Rt))
    raw.Ht<-matrix(ncol=length(object$ord.time),nrow=nrow(object$rep.Rt))

    for(i in 1:nrow(object$rep.Rt)){
	finite<-is.finite(object$rep.Rt[i,])
      low.lm<-min(object$rep.Rt[i,finite])-100
      Rt.fun<-stepfun(object$ord.time[finite],c(low.lm,object$rep.Rt[i,finite]))
      new.Ht[i,]<-Rt.fun(arg$new.time)
      raw.Ht[i,]<-Rt.fun(object$ord.time)
    }
    if(p==0) rep.bz<-0
    if(p>0) rep.bz<-object$rep.beta%*%arg$newdata
    new.error<-new.Ht+matrix(rep.bz,ncol=n,nrow=nrow(object$rep.Rt))
    e.sd<-apply(new.error,2,sd)
    ul.err<-new.err+qnorm(1-arg$alpha/2)*e.sd
    ll.err<-new.err-qnorm(1-arg$alpha/2)*e.sd
    ul.err<-monot(ul.err)
    ll.err<-monot(ll.err)

    raw.error<-raw.Ht+matrix(rep.bz,ncol=length(object$ord.time),nrow=nrow(object$rep.Rt))
    Qa<-quantile(apply(t(na.omit(t(abs(scale(raw.error))))),1,max),1-arg$alpha,na.rm=TRUE)

    ub.err<-new.err+Qa*e.sd
    lb.err<-new.err-Qa*e.sd
    ub.err<-monot(ub.err)
    lb.err<-monot(lb.err)

    if(r==0){
      ll.st<-exp(-exp(ul.err))
      ul.st<-exp(-exp(ll.err))

      lb.st<-exp(-exp(ub.err))
      ub.st<-exp(-exp(lb.err))
    } else{
      ll.st<-exp(-log(1+r*exp(ul.err))/r)
      ul.st<-exp(-log(1+r*exp(ll.err))/r)

      lb.st<-exp(-log(1+r*exp(ub.err))/r)
      ub.st<-exp(-log(1+r*exp(lb.err))/r)
    }

    ll.st[ll.st<0]<-lb.st[lb.st<0]<-0
    ul.st[ul.st>1]<-ub.st[ub.st>1]<-1
    pred<-data.frame(pred,ll.st,ul.st,lb.st,ub.st)
    pred$Qa=Qa
    }

  class(object)<-"TransModel"
  class(pred)<-"predict.TransModel"
  return(pred)
}
