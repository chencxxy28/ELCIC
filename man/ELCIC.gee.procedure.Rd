\name{ELCIC.gee.procedure}
\alias{ELCIC.gee.procedure}
\title{
The whole procedure for joint selection of mean structure and correlation structure in longitudinal data
}
\description{
This function provides the overall procedure for joint selection of mean structure and correlation structure in longitudinal data
}
\usage{
ELCIC.gee.procedure(x,y,r,id,time,candidate.sets,name.var=NULL,dist,candidate.cor.sets=c("independence","exchangeable", "ar1"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{r}{a vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.}
  \item{id}{a vector indicating subject id.}
  \item{time}{the number of observations for each subject}
  \item{candidate.sets}{a list containing index corresponding to candidate covariates}
  \item{name.var}{a list containing names of candidate covariates. The names should be subset of column names of x matrix.}
  \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
  \item{candidate.cor.sets}{a vector containing condidate correlation structures. It is c("independence","exchangeable", "ar1") as default}
}
\value{
A matrix with each element containing ELCIC value for each candidate model
}
