S0.fun <-
function(ord.bz,n,ord.wt){
  s0 = numeric(n)
  a = ord.wt*exp(-ord.bz)
  s0 = sum(a) - cumsum(a) + a
  return(s0)
}
