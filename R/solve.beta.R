solve.beta <-
function(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max){
    p = length(beta.ini)
    n = length(Rt.ini)

    if(p>1) zc = t(t(z) - apply(z,2,mean))
    if(p==1) zc = as.matrix(z-mean(z))
    ix = order(obs.t)
    ord.delta = delta[ix]
    ord.z = as.matrix(zc[ix,])
    colnames(ord.z)= colnames(z)
    ord.wt = wt[ix]
    
    ###solving estimating equations for beta and Rt
    dif = 1
    iter = 0
    while(dif > dx & iter <= iter.max){
        ord.bz = as.numeric(ord.z%*%beta.ini)
        if(r == 0){
           s0 = S0.fun(ord.bz,n,ord.wt)
           s1 = S1.fun(ord.z,ord.bz,n,p,ord.wt)
           Ubeta = Ubeta1.fun(ord.delta,ord.z,s0,s1,ord.wt)
           Hbeta = Hbeta1.fun(ord.delta,ord.z,ord.bz,s0,s1,p,ord.wt)
	     Rt = Rt.fun(ord.delta,s0,ord.wt)
	     Rt.ini = Rt
        }
        if(r > 0){
           s0 = S02.fun(ord.bz,Rt.ini,r,n,ord.wt)
           Rt = Rt.fun(ord.delta,s0,ord.wt)
           Ubeta = Ubeta2.fun(ord.delta,ord.z,ord.bz,Rt,r,ord.wt)
           Hbeta = Hbeta2.fun(ord.z,ord.bz,Rt,r,p,ord.wt)
           Rt.ini = Rt
        }
        beta = beta.ini + as.numeric(ginv(Hbeta)%*%Ubeta)  
        
        dif = max(abs(beta-beta.ini))   
        iter = iter + 1
        beta.ini = beta
    }
    converged = as.numeric(dif > dx)   
    return(list(Beta=beta,converged=converged,iter=iter,data=data.frame(ord.time=sort(obs.t),ord.delta=ord.delta,Rt=Rt),Hbeta=Hbeta))
}
