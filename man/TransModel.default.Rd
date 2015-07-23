\name{TransModel.default}
\alias{TransModel.default}

\title{
Default method for function TransModel.
}
\description{
Default method for function TransModel.
}
\usage{
\method{TransModel}{default}(formula = formula(data), data = parent.frame(),r=0,
CICB.st=FALSE,dx=0.001,iter.max=100,num.sim=200)
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
