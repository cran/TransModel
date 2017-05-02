Ubeta1.fun <-
function(ord.delta,ord.z,S0,S1,ord.wt){
  a = -ord.delta*(ord.z-S1/S0)*ord.wt
  b = apply(a,2,sum)
  return(b)  
}
