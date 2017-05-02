monot <-
function(aa){ while(sum(diff(c(-Inf,aa))<0)>0) aa[diff(c(-Inf,aa))<0]<-aa[which(diff(c(0,aa))<0)-1];return(aa)}
