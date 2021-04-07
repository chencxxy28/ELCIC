\name{lambda.find.gee}
\alias{lambda.find.gee}
\title{
calculate the tuning parameters involved in ELCIC
}
\description{
This function provides an efficient algorithm to calculate the tuning parameters involved in ELCIC.
}
\usage{
lambda.find.gee(x, y, id, betahat, r, dist, ro, phi, corstr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{r}{a vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.}
  \item{id}{a vector indicating subject id.}
  \item{betahat}{a plug-in estimator solved by an external estimation procedure, such as GEE.}
  \item{ro}{a correlation coefficients obtained from an external estimation procedure, such as GEE.}
  \item{phi}{an over-dispersion parameter obtained from an external estimation procedure, such asGEE}
   \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
   \item{corstr}{a condidate correlation structure. It can be "independence","exchangeable", and "ar1"}

}
\value{
tuning parameter values
}
