var.e <-
function(obs.t,delta,z,r,dx,iter.max,n,p,num.sim){
  beta.est = matrix(nrow=num.sim,ncol=p)
  Rt.est <-Rt0.est <- matrix(nrow=num.sim,ncol=n)
  i<-1
  while(i<=num.sim){
	wt = rexp(n)
	beta.ini = rep(0,p)
	Rt.ini = cumsum(1/(n:1))
	temp = try(solve.beta(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max),silent=TRUE)
      if(!is.character(temp)){
	 succ1<-sum(is.na(c(temp$Beta,temp$data$Rt)))+sum(!is.finite(c(temp$Beta,temp$data$Rt)))
	 if(succ1==0){
	  beta.est[i,] = temp$Beta
	  Rt.est[i,] = temp$data$Rt
        if(p>1) Rt0.est[i,]<-log(temp$data$Rt)-sum(temp$Beta*apply(z,2,mean))
        if(p==1) Rt0.est[i,]<-log(temp$data$Rt)-sum(temp$Beta*mean(z))
        i<-i+1
	 }
      }
  }
return(list(beta.est,Rt0.est))
}
