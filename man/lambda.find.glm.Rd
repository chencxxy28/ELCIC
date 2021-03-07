\name{lambda.find.glm}
\alias{lambda.find.glm}
\title{
The tuning parameter in empirical likelihood
}
\description{
This function aims to efficiently calculate the tuning parameter lambda in the empirical likelihood.
}
\usage{
lambda.find.glm(x, y, betahat, dist)
}
\arguments{
  \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{betahat}{a plug-in estimator solved by an external estimating procedure.}
 \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
}
\value{
A value of lambda involved in the empirical likelihood.
}
