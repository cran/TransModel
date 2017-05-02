\name{plot.predict.TransModel}
\alias{plot.predict.TransModel}

\title{
Plot the predicted survival curve with confidence interval/band.
}
\description{
Survival curve with confidence interval or confidence band will be plotted for the returned object from predict.TransModel.
}
\usage{
\method{plot}{predict.TransModel}(x, CI = FALSE, CB = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object from the function predict.TransModel.
}
  \item{CI}{
Whether a 100(1-alpha)\% pointwise confidence interval will be added to the survival curve.
}
  \item{CB}{
Whether a 100(1-alpha)\% confidence band will be added to the survival curve.
}
  \item{\dots}{
Other plot arguments.
}
}
\note{
Either CI or CB is set to be TRUE in this plot method, the object from predict.TransModel needs to have the
 argument CICB.st to be TRUE as well. Otherwise, no confidence interval or band will be added to the plot.
}


