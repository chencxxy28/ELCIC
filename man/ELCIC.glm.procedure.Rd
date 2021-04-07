\name{ELCIC.glm.procedure}
\alias{ELCIC.glm.procedure}
\title{The whole variable selection procedure in GLM}
\description{
This function provides the overall procedure for variable selection in GLM
}
\usage{
ELCIC.glm.procedure(x,y,candidate.sets,name.var=NULL,dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{candidate.sets}{a list containing index corresponding to candidate covariates}
  \item{name.var}{a list containing names of candidate covariates. The names should be subset of column names of x matrix.}
  \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
}
\value{
A matrix with each element containing ELCIC value for each candidate model (in columns) and (in rows)

}
