\name{gee.generator}
\alias{gee.generator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate longitudinal data without missingness}
\description{
A function for generating longitudinal data without missingness
}
\usage{
gee.generator(beta,samplesize,time,num.time.dep,num.time.indep,rho=0.5,x.rho=0.2,dist,cor.str="exchangeable",x.cor.str="exchangeable")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{a vector containing underlying true coefficients for each covariate in the model (including the intercept).}
  \item{samplesize}{the sample size.}
  \item{time}{the number of observations per subject.}
  \item{num.time.dep}{the number of time-dependent covariates.}
  \item{num.time.indep}{the number of time-independent covariates.}
  \item{rho}{the correlation coefficient for residuals across time.}
  \item{x.rho}{the correlation coefficient for time-dependent covariates across time.}
 \item{dist}{a specified distribution. It can be "gaussian", "poisson",and "binomial".}
 \item{cor.str}{the correlation structure for residuals across time. It can be "independence","exchangeable", and "ar1".}
 \item{x.cor.str}{the correlation structure for time-dependent covariates across time. It can be "independence","exchangeable", and "ar1".}
}

\value{
\item{x}{a matrix containing covariates. The first column should contain all ones corresponding to the intercept.}
  \item{y}{a vector containing outcomes.}
  \item{id}{a vector indicating subject id.}
}
