\name{ELCIC.glm}
\alias{ELCIC.glm}
\title{
Variable selection in generalized linear models
}
\description{
This function will provide values of several model selection criteria including AIC, BIC, GIC, and ELCIC.
}
\usage{
ELCIC.glm(x, y, index.var, name.var = NULL, dist)
}
\arguments{
  \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{index.var}{a vector containing index corresponding to candidate covariates.}
  \item{name.var}{a vector containing names of candidate covariates. The names should be subset of column names of x matrix.}
  \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
}
\value{
\item{ELCIC}{a value of ELCIC.}
\item{AIC}{a value of AIC}
\item{BIC}{a value of BIC}
\item{GIC}{a value of GIC}
}

