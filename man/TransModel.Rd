\name{TransModel}
\alias{TransModel}

\title{
Fit linear transformation models for censored data
}
\description{
This function is used to fit a linear transformation model, such as the proportional hazards model or 
proportional odds model, to the right censored survival data.
}
\usage{
TransModel(formula = formula(data), data = parent.frame(),r,CICB.st,dx,iter.max,num.sim)
}
\arguments{
  \item{formula}{
A survival formula based on Surv function, containg survival time, right censoring indicator and coavariates.
}
  \item{data}{
Data set with all the variables needed in formula.
}
  \item{r}{
Parameter in the hazard function, used to define different linear models. See details for more information.
}
  \item{CICB.st}{
Whether or not the perturbation for the 95\% confidence interval and confidence bands will be done. The default value is FALSE.  
}
  \item{dx}{
Upbound of convergence limit. Default is 0.001.
}
  \item{iter.max}{
Maximum number of iterations before convergence. Default is 100.
}
  \item{num.sim}{
The number of perturbation, only works when CICB.st=TRUE. Default is 200.
}
}
\details{
In the linear transformation model H(t)=-b'z+e, the hazard function for error term e is defined as:
h(x)=exp(x)/(1+r*exp(x)),
where the parameter r can be changed for different models. For example, r=0 refers to the proportional hazards
model and r=1 refers to a proportional odds model. The default value for r is 0.
}
\value{
\item{coefficients}{Estimated coefficients for covariates in the specified linear transformation model.}
\item{vcov}{Estimated covariance matix for the coefficients.}
\item{ord.time}{Sorted survival time of original dataset.}
\item{St}{Baseline survival probabilities with regard to ord.time.}
\item{converged}{Convergence status, 0 indicates converged, and number of iterations used for convergence.}
}
\references{
Kani Chen, et al., Semiparametric analysis of transformation models with censored data. Biometrika, 
89(3), 659-668, 2002.
}
\examples{
beta0 = c(1,-1)
# Fit proportional hazards model
data(PH_examp)
mod1<-TransModel(formula=Surv(time,status)~gender+age,data=PH_examp,r=0)
print(mod1)
summary(mod1)
mod1$coefficients
mod1$vcov
mod1$converged

# Estimated baseline survival function
pred1<-predict(mod1,covar=data.frame(gender=0,age=0))
plot(pred1)
lines(pred1$time,exp(-((exp(2*pred1$time)-1)/2)),col=3) # True value

### Not Run ###
# survival estimate with 95% pointwise CI and overall CB
# mod1<-TransModel(formula=Surv(time,status)~gender+age,data=PH_examp,r=0,CICB.st=TRUE,num.sim=50)
# pred1<-predict(mod1,covar=data.frame(gender=0,age=0))
# plot(pred1,lty=1,col=1,CI=TRUE,CB=TRUE)
# lines(pred1$time,exp(-((exp(2*pred1$time)-1)/2)),col=3) # True value

# compare survival curve for gender=0 or 1, age=0
# pred2<-predict(mod1,covar=data.frame(gender=1,age=0))
# plot(pred1,lty=2,col=1)
# lines(pred2$time,pred2$survival,type="s",lty=2,col=2)
# lines(pred1$time,exp(-((exp(2*pred1$time)-1)/2)),col=1) # True value
# lines(pred2$time,exp(-((exp(2*pred2$time)-1)/2))^exp(sum(beta0*c(1,0))),col=2) # True value

# Fit proportional odds model
data(PO_examp)
mod2=TransModel(Surv(time,status)~gender+age,data=PO_examp,r=1)
print(mod2)
summary(mod2)
pred2<-predict(mod2,covar=data.frame(gender=0,age=0))
}
\keyword{ Linear transformation models }
