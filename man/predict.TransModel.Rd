\name{predict.TransModel}
\alias{predict.TransModel}
\title{
Predict survival probabilities for the data or at 
new given time points.
}
\description{
Calculate survival probabilities, 95 percent confidence intervals and confidence bands from the estimated 
transformation model. If new time points is not specified survival probabilities of the original data set 
will be returned. If covariate variables are not specified the baseline survival probabilities will be calculated.}
\usage{
\method{predict}{TransModel}(object, ...)
}
\arguments{
  \item{object}{
An object from the function TransModel.
}
  \item{...}{
See details.
}
}

\details{ Other arguments specified in ... including:
covar: data frame of the values for each covariate variables. If not specified, 0 will be used for all variables.
new.time: Vector of ordered time points to be used for survival probability calculation. If null, original time points
in the data set will be used.
}
\value{
Data frame with time points and predicted survival probabilities will be returned. If CICB.st is TRUE, the lower
 limits and upper limits of the confidence interval and confidence bands are also included in the data frame.
}

