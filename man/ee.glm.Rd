\name{ee.glm}
\alias{ee.glm}
\title{
Estimating equations for ELCIC under GLM
}
\description{
A specified estimating equations for ELCIC under GLM
}
\usage{
ee.glm(x, y, betahat, dist)
}
\arguments{
  \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{betahat}{a plug-in estimator solved by an external estimating procedure.}
 \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
}
\value{
A matrix containing values of calculated estimating equations.
}
