\name{summary.TransModel}
\alias{summary.TransModel}
\title{
Summary results for fitting the linear transformation model
}
\description{
A summary table contains the coefficient estimates, standard errors, test statistics and p values in the linear
 transformation model. The loglikelihood value and the AIC value are also presented at the bottom of the table.
}
\usage{
\method{summary}{TransModel}(object, ...)
}

\arguments{
  \item{object}{
An object returned from the function TransModel.
}
  \item{\dots}{
Other auguments to be specified for the summary function.
}
}

\value{
\item{coefficients}{The summary table.}
\item{loglik}{The loglikelihood value of the fitted model.}
\item{AIC}{The AIC value.}
}

\examples{
data(veteran)
fit<-TransModel(Surv(time,status)~karno+as.factor(celltype),data=veteran,r=0,subset=(prior==0))
summary(fit)

fit0<-TransModel(Surv(time,status)~1,data=veteran,r=0,subset=(prior==0))  #the null model
summary(fit0)

#likelihood ratio test
LRT<--2*(summary(fit0)$loglik-summary(fit)$loglik)
df<-fit$p-fit0$p
1-pchisq(LRT,df)
}

