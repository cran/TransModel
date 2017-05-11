\name{summary.TransModel}
\alias{summary.TransModel}
\title{
Summary results for fitting the linear transformation model
}
\description{
A summary table contains the coefficient estimates, standard errors, test statistics and p values in the linear
 transformation model.
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
}

\examples{
data(veteran)
fit<-TransModel(Surv(time,status)~karno+as.factor(celltype),data=veteran,r=0,subset=(prior==0))
summary(fit)

fit0<-TransModel(Surv(time,status)~1,data=veteran,r=0,subset=(prior==0))  #the null model
summary(fit0)
}

