\name{glm.generator}
\alias{glm.generator}
\title{
Data generation under GLMs
}
\description{
A function provides simulated outcomes as well as responses under the framework of GLMs.
}
\usage{
glm.generator(beta, samplesize, rho = 0, dist, sd.gaussian = NULL, ov = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{The underlying true coefficient for each covariates in the model (including the intercept).}
  \item{samplesize}{the sample size.}
  \item{rho}{the correlation coefficient among covariates.}
 \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
  \item{sd.gaussian}{the standard deviation for the outcome from Gaussian distribution.}
  \item{ov}{the dispersion parameter for the outcome from Negative Binomial distribution.}
}
\value{
  \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
}
