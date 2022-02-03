plot.predict.TransModel <-
function(x,CI=FALSE,CB=FALSE,...){
  arg<-list(...)
  if(is.null(arg$xlab)) arg$xlab<-"Time"
  if(is.null(arg$ylab)) arg$ylab<-"Survival Probability"
  if(is.null(arg$ylim))  arg$ylim<-c(0,1)
  plot(x$time,x$survival,type="s",xlab=arg$xlab,ylab=arg$ylab,...)
  if(CI){
     lines(x$time,x$low.ci,lty=2,type="s")
     lines(x$time,x$up.ci,lty=2,type="s")
  }
  if(CB){
     lines(x$time,x$low.cb,lty=3,type="s")
     lines(x$time,x$up.cb,lty=3,type="s")
  }
  class(x)<-"predict.TransModel"
}
