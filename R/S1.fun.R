S1.fun <-
function(ord.z,ord.bz,n,p,ord.wt){
  s1 = matrix(0,nrow=n,ncol=p)
  a = ord.z*ord.wt*exp(-ord.bz)
  for(i in 1:p){
    temp = as.numeric(a[,i])
    s1[,i] = sum(temp) - cumsum(temp) + temp
  }
  return(s1)
}
