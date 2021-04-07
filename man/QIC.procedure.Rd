\name{QIC.procedure}
\alias{QIC.procedure}
\title{
Joint selection of marginal mean and correlation structures in longitudinal data based on QIC
}
\description{
This function provides the Joint selection of marginal mean and correlation structures in longitudinal data based on QIC
}
\usage{
QIC.procedure(x,y,id,dist,candidate.sets,name.var=NULL,candidate.cor.sets)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{id}{a vector indicating subject id.}
  \item{candidate.sets}{a list containing index corresponding to candidate covariates}
  \item{name.var}{a list containing names of candidate covariates. The names should be subset of column names of x matrix.}
  \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
  \item{candidate.cor.sets}{a vector containing condidate correlation structures. It is c("independence","exchangeable", "ar1") as default}
}
\value{
A vector with each element containing QIC value for each candidate model. The row name of this vector is the selected correlation structure.
}
