\name{ELCIC.gee}
\alias{ELCIC.gee}
\title{
calculate ELCIC value for a given candidate model
}
\description{
This is the function to calculate ELCIC value for a given candidate model
}
\usage{
ELCIC.gee(x, y, r, id, time, index.var, name.var = NULL, dist, corstr)
}
\arguments{
 \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{r}{a vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.}
  \item{id}{a vector indicating subject id.}
  \item{time}{the number of observations for each subject}
  \item{index.var}{a vector containing index corresponding to candidate covariates.}
  \item{name.var}{a vector containing names of candidate covariates. The names should be subset of column names of x matrix.}
  \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
  \item{corstr}{a candidate correlation structure. It can be "independence","exchangeable", and "ar1".}
}
\value{
A ELCIC value for a given candidate model.
}
