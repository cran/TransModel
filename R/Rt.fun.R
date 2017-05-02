Rt.fun <-
function(ord.delta,s0,ord.wt){
    Rt = cumsum(ord.wt*ord.delta/s0)
    return(Rt)
}
