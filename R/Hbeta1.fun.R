Hbeta1.fun <-
function(ord.delta,ord.z,ord.bz,S0,S1,p,ord.wt){
  a = matrix(0,nrow=p,ncol=p)
  for(j in 1:p){
    for(k in 1:p){
      temp = as.numeric(ord.z[,j])*as.numeric(ord.z[,k])*exp(-ord.bz)*ord.wt
      b = sum(temp) - cumsum(temp) + temp
      a[j,k] = sum(ord.wt*ord.delta*b/S0)
    }
  }
  e1 = ord.delta*S1/S0*ord.wt
  e2 = S1/S0
  a = a - t(e1)%*%e2
  return(a)
}
