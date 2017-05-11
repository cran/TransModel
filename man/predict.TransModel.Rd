\name{predict.TransModel}
\alias{predict.TransModel}

\title{
Predict survival probabilities with confidence intervals and confidence bands
}
\description{
Calculate survival probabilities, confidence intervals and confidence bands in the estimated 
transformation model. If new time points is not specified, event times in the original data set 
will be used. If no new covariate values are specified, zeros will be used and the baseline 
survival probabilities will be calculated.
}
\usage{
\method{predict}{TransModel}(object, ...)
}

\arguments{
  \item{object}{
An object returned from the function TransModel.
}
  \item{\dots}{
Other arguments including:

newdata: a vector containing the values for each covariate variables specified in the model. If not specified, 0 will be used for all variables.

new.time: vector of ordered time points to be used for survival probability calculation. If null, distinct event time points
in the original dataset will be used.

alpha: used to determine the confidence level of the predicted confidence interval/band for the survival curve. 
The default value is 0.05, corresponding to a 95\% confidence level.
}
}

\value{
\item{time}{ordered time points on which survival probabilities are calculated.}
\item{survival}{predicted survival probabilities.}
\item{low.ci}{the lower limit of confidence interval.}
\item{up.ci}{the upper limit of confidence interval.}
\item{low.cb}{the lower limit of confidence band. }
\item{up.cb}{the upper limit of confidence band.}
}

\note{
The values low.ci, up.ci, low.cb and up.cb will be returned only if CICB.st=TRUE in the original model.
}


\examples{
data(PH_examp)
mod1<-TransModel(formula=Surv(time,status)~gender+age,data=PH_examp,r=0)
# Predict survival probability
pred1<-predict(mod1,newdata=c(0,1))
pred2<-predict(mod1,newdata=c(1,1))
plot(pred1)
lines(pred2$time,pred2$survival,type="s",col=2)

### Not Run ###
# survival estimate with 95% pointwise CI and overall CB
# mod1<-TransModel(formula=Surv(time,status)~gender+age,data=PH_examp,r=0,CICB.st=TRUE,num.sim=50)
# pred1<-predict(mod1,newdata=c(0,1))
# plot(pred1,lty=1,col=1,CI=TRUE,CB=TRUE)
# Change the confidence level to 90%
# pred1<-predict(mod1,newdata=c(0,1),alpha=0.1)
}

