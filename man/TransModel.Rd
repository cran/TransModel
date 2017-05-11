\name{TransModel}
\alias{TransModel}
\title{
Main function for fitting the linear transformation models for right censored data
}
\description{
This function is used to fit a linear transformation model, such as the proportional hazards model or 
proportional odds model, to the right censored survival data.
}
\usage{
TransModel(formula = formula(data), data = parent.frame(), r, CICB.st, subset, 
dx, iter.max, num.sim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A survival formula based on the Surv() function, containg survival time, right censoring indicator and coavariates.
}
  \item{data}{
Data set with all the variables needed in formula.
}
  \item{r}{
Parameter in the hazard function, used to define different linear models. See details for more information.
}
  \item{CICB.st}{
Whether or not the perturbation for deriving the confidence intervals and confidence bands of survival estimates will be done. The default value is FALSE.  
}
  \item{subset}{
Conditions for subsetting the dataset.
}
  \item{dx}{
Convergence tolerance. Default is 0.001.
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
h(x)=exp(x)/(1+r*exp(x)), where the parameter r must be a non-negative value and can be changed for different models. For example, r=0 refers to the proportional hazards
model and r=1 refers to a proportional odds model. The default value for r is 0.
}
\value{
\item{coefficients}{Estimated coefficients for covariates in the specified linear transformation model.}
\item{vcov}{Estimated covariance matix for the coefficients.}
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

# Fit proportional odds model
data(PO_examp)
mod2=TransModel(Surv(time,status)~gender+age,data=PO_examp,r=1)
print(mod2)
summary(mod2)
}
\keyword{ Linear transformation models }
