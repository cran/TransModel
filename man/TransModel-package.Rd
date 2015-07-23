\name{TransModel-package}
\alias{TransModel-package}
\docType{package}
\title{
Fit linear transformation models for censored data
}
\description{
Semiparametric analysis for linear transformation model, such as the proportional hazards model
 and the proportional odds model, is fitted using Chen's method for right censored survival data.
}
\details{
\tabular{ll}{
Package: \tab TransModel\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-23\cr
License: \tab GPL-2\cr
Depends: stats,graphics,survival
}
The main function TransModel gives estimates for coefficients and covariance matrix for the predictors
in the linear transformation model. Print and summary method can be applied to the returned object.
}
\author{
Jie Zhou, Jiajia Zhang, Wenbin Lu
Maintainer: Jie Zhou <zhoujie02569@gmail.com>
}
\references{
Kani Chen, et al., Semiparametric analysis of transformation models with censored data. Biometrika, 
89(3), 659-668, 2002.
}
