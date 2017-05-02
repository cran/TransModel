plot.predict.TransModel <-
function(x,CI=FALSE,CB=FALSE,...){
  arg<-list(...)
  if(is.null(arg$xlab)) arg$xlab<-"Time"
  if(is.null(arg$ylab)) arg$ylab<-"Survival Probability"
  if(is.null(arg$ylim))  arg$ylim<-c(0,1)
  plot(x$time,x$survival,type="s",xlab=arg$xlab,ylab=arg$ylab,...)
  if(CI){
     lines(x$time,x$ll.st,lty=2,type="s")
     lines(x$time,x$ul.st,lty=2,type="s")
  }
  if(CB){
     lines(x$time,x$lb.st,lty=3,type="s")
     lines(x$time,x$ub.st,lty=3,type="s")
  }
  class(x)<-"predict.TransModel"
}
