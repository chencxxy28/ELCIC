##calculate I  #Now works for gaussian, poisson, and binomial distribution
II<-function (x,y,beta,size,samplesize,dist)
{
    I_G<-matrix(0,nrow=size,ncol=size)
    if(dist=="gaussian")
    {
    for (i in 1:samplesize)
    {
        fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%beta))
        I_G<-I_G+fitted%*%t(fitted)
    }
    }else if(dist=="poisson")
    {
        for (i in 1:samplesize)
        {
            fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%beta))
            I_G<-I_G+fitted%*%t(fitted)
        }
    }else if(dist=="binomial")
    {
        for (i in 1:samplesize)
        {
            fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-1/(exp(-t(as.matrix(x[i,],ncol=1))%*%beta)+1))
            I_G<-I_G+fitted%*%t(fitted)
        }
    }
    I_G
}

#calcuate J #Now works for gaussian, poisson, and binomial distribution
JJ<-function (x,y,beta,size,samplesize,dist)
{
    J_G<-matrix(0,nrow=size,ncol=size)
    if(dist=="gaussian")
    {
    for (i in 1:samplesize)
    {
        fitted<-as.matrix(x[i,],ncol=1)%*%t(as.matrix(x[i,],ncol=1))
        J_G<-J_G+fitted
    }
    }else if(dist=="poisson")
    {
        for (i in 1:samplesize)
        {
            fitted<-as.vector(exp(t(as.matrix(x[i,],ncol=1))%*%beta))*as.matrix(x[i,],ncol=1)%*%t(as.matrix(x[i,],ncol=1))
            J_G<-J_G+fitted
        }
    }else if (dist=="binomial")
    {
        for (i in 1:samplesize)
        {
            pi_i<-1/(exp(-t(as.matrix(x[i,],ncol=1))%*%beta)+1)
            fitted<-as.vector(pi*(1-pi))*as.matrix(x[i,],ncol=1)%*%t(as.matrix(x[i,],ncol=1))
            J_G<-J_G+fitted
        }
    }
    J_G
}

#'@title Cross-sectional data generation under GLM
#'@description A function provides simulated outcomes as well as covariates under the framework of GLM. All covariates (except for intercept) are normally distributed.
#'@usage glm.generator(beta, samplesize, rho = 0, dist, sd.gaussian = NULL, ov = NULL)
#'@param beta The underlying true coefficient for each covariates in the model (including the intercept).
#'@param samplesize The sample size.
#'@param rho The correlation coefficient among covariates.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param sd.gaussian The standard deviation for the outcome from Gaussian distribution.
#'@param ov The dispersion parameter for the outcome from Negative Binomial distribution.
#'
#'
#'@return x: a matrix containing continuous covariates. The first column should contain all ones corresponding to the intercept.
#'@return y: a vector containing outcomes.
#'
#'@examples
#'beta<-c(0.5,0.5,0.5,0)
#'samplesize<-100
#'data<-glm.generator(beta=beta,samplesize=samplesize,rho=0.5,dist="poisson")
#'
#'@export
#'@importFrom mvtnorm rmvnorm
#'@import MASS

#generate glm data
glm.generator<-function(beta,samplesize,rho=0,dist,sd.gaussian=NULL,ov=NULL)
{
    pn<-length(beta)-1
    X.sigma <- matrix(0,(pn),(pn)) #covariance for the covariate matrix
    {
        for (i in 1:(pn))
            X.sigma[i,] <- rho^(abs((1:(pn))-i))
    }
    xtrue<-rmvnorm(samplesize, mean = rep(0,pn), X.sigma) #generate a covariate matrix
    x<-cbind(rep(1,samplesize),xtrue) #the overall covariate matrix with intercept

    switch(dist,
           "gaussian"={y <- rnorm(samplesize,(x%*%beta),sd.gaussian)
           },
           "binomial"={y <- rbinom(n=samplesize,size=1,prob=1/(1+exp(-x%*%beta)))
           },
           "poisson"={y <- rpois(samplesize,exp(x%*%beta))
           },
           "NB"={y <- rnbinom(n=samplesize,size=ov,mu=exp(x%*%beta))
           },
           stop("Invalid type of dist. It should be one of gaussian,binomial,poisson,NB")
    )

    # if(dist=="gaussian")
    # {
    # y<-rnorm(samplesize,(x%*%beta),sd.gaussian)
    # }else if(dist=="poisson")
    # {
    #     y<-rpois(samplesize,exp(x%*%beta))
    # }else if(dist=="binomial")
    # {
    #     y<-rbinom(n=samplesize,size=1,prob=1/(1+exp(-x%*%beta)))
    # }else if(dist=="NB")
    # {
    #     y<-rnbinom(n=samplesize,size=ov,mu=exp(x%*%beta))
    # }
    return(list(y=y,x=x))
}

#'@title Generate longitudinal data without missingness
#'@description A function for generating longitudinal data without missingness. All covariates (except for intercept) are normally distributed.
#'@usage gee.generator(beta,samplesize,time,num.time.dep,num.time.indep,
#'      rho,x.rho,dist,cor.str,x.cor.str)
#'@param beta A vector containing underlying true coefficients for each covariate in the model (including the intercept).
#'@param samplesize The sample size.
#'@param time The number of observations per subject.
#'@param num.time.dep The number of time-dependent covariates.
#'@param num.time.indep The number of time-independent covariates (not include intercept).
#'@param rho The correlation coefficient for residuals across time.
#'@param x.rho The correlation coefficient for time-dependent covariates across time.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param cor.str The correlation structure for residuals across time. It can be "independence","exchangeable", and "ar1".
#'@param x.cor.str The correlation structure for time-dependent covariates across time. It can be "independence","exchangeable", and "ar1".
#'
#'@return x: a matrix containing continuous covariates. The first column should contain all ones corresponding to the intercept.
#'@return y: a vector containing outcomes.
#'@return id: a vector indicating subject id.
#'
#'@examples
#'beta<-c(-1,1,0.5,0)
#'samplesize<-100
#'geesimdata<-gee.generator(beta=beta,samplesize=samplesize,time=3,num.time.dep=2,
#'num.time.indep=1,rho=0.4,x.rho=0.2,dist="poisson",cor.str="exchangeable",
#'x.cor.str="exchangeable")
#'geesimdata$y
#'
#'@export
#'@importFrom mvtnorm rmvnorm
#'@importFrom PoisNor genPoisNor
#'@importFrom bindata rmvbin
#'
gee.generator<-function(beta,samplesize,time,num.time.dep,num.time.indep,rho,x.rho,dist,cor.str,x.cor.str) {

    if(length(beta)!=(num.time.dep+num.time.indep+1))
        {
        stop("Invalid input for number of time independent covariate and dependent covariate. Sum should be equal to the length of beta minus one")
        }
    switch(cor.str,
           "exchangeable"={
               sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
               {
                   for (i in 1:(time))
                       sigma.dep[i,] <- rho
               }
           },
           "ar1"={
               sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
               {
                   for (i in 1:(time))
                       sigma.dep[i,] <- rho^(abs((1:(time))-i))
               }
           },
           "independence"={
               sigma.dep <- matrix(0,(time),(time))
           },
           stop("Invalid type of correlation.struture for outcomes. It should be one of exchangeable, ar1, and independence")
    )
    diag(sigma.dep)<-1

    # if(cor.str=="exchangeable")
    # {
    #     sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    #     {
    #         for (i in 1:(time))
    #             sigma.dep[i,] <- rho
    #     }
    # }else if(cor.str=="ar1"){
    #     sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    #     {
    #         for (i in 1:(time))
    #             sigma.dep[i,] <- rho^(abs((1:(time))-i))
    #     }
    # }else{
    #     sigma.dep <- matrix(0,(time),(time))
    # }
    # diag(sigma.dep)<-1

    switch(x.cor.str,
           "exchangeable"={
               x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
               {
                   for (i in 1:(time))
                       x.sigma.dep[i,] <- x.rho
               }
           },
           "ar1"={
               x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
               {
                   for (i in 1:(time))
                       x.sigma.dep[i,] <- x.rho^(abs((1:(time))-i))
               }
           },
           "independence"={
               x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
           },
           stop("Invalid type of correlation.struture for covariates. It should be one of exchangeable, ar1, and independence")
    )
    diag(x.sigma.dep)<-1

    # if(x.cor.str=="exchangeable")
    # {
    #     x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    #     {
    #         for (i in 1:(time))
    #             x.sigma.dep[i,] <- x.rho
    #     }
    # }else if(x.cor.str=="ar1"){
    #     x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    #     {
    #         for (i in 1:(time))
    #             x.sigma.dep[i,] <- x.rho^(abs((1:(time))-i))
    #     }
    # }else{
    #     x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    # }
    # diag(x.sigma.dep)<-1

    xf<-rep()
    y<-rep()

    #full X
    for (jj in 1:samplesize)
    {
    x.dep<-t(rmvnorm(num.time.dep,rep(0,time=time),x.sigma.dep))
    x.indep<-rep(rnorm(num.time.indep,0,1),time=time)
    xf.i<-cbind(1,x.dep,x.indep)
    #true X
    switch(dist,
           "gaussian"={
               mu<-xf.i%*%beta
               y.i<-rmvnorm(1,mu,sigma.dep)
           },
           "binomial"={
               mu<-1/(1+exp(-xf.i%*%beta))
               y.i<-rmvbin(n=1,margprob=mu,bincorr=sigma.dep)
           },
           "poisson"={
               mu<-exp(xf.i%*%beta)
               y.i<-genPoisNor(n=1,no.pois=length(mu),lamvec=mu,cmat.star=sigma.dep,no.norm=0,mean.vec=NULL,sd.vec=NULL)
           },
           stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")
    )
    # if(dist=="gaussian")
    # {
    #     mu<-xf.i%*%beta
    #     y.i<-rmvnorm(1,mu,sigma.dep)
    # }else if(dist=="poisson"){
    #     mu=exp(xf.i%*%beta)
    #     y.i<-genPoisNor(n=1,no.pois=length(mu),lamvec=mu,cmat.star=sigma.dep,no.norm=0,mean.vec=NULL,sd.vec=NULL)
    # }else if(dist=="binomial"){
    #     mu=1/(1+exp(-xf.i%*%beta))
    #     y.i<-rmvbin(n=1,margprob=mu,bincorr=sigma.dep)
    # }
    xf<-rbind(xf,xf.i)
    y<-c(y,y.i)
    }
    id<-rep(1:samplesize,each=time)
    return(list(y=y,x=xf,id=id))
}



v<-function(x,beta,dist) #function used to calculate variance, mean, and derivative of random variables followed by binomial, poisson, and gaussian
{
    n<-length(x[,1])
    reg<-x%*%beta
    pi<-exp(reg)/(1+exp(reg))
    if (dist=="gaussian")
    {
        v<-rep(1,n)
        der<-x
        mu<-reg
    }
    else if (dist=="binomial")
    {
        v<-pi*(1-pi)
        der<-x*as.vector(v)
        mu<-as.vector(pi)
    }
    else if (dist=="poisson")
    {
        v<-exp(reg)
        der<-x*as.vector(exp(reg))
        mu<-as.vector(exp(reg))
    }
    list(v=v,der=der,mu=mu)
}


R<-function(rho,id) #working correlation matrix
{
    time<-length(id)/length(unique(id))
    R<-diag(1,time,time)
    rho<-c(1,rho)
    for (i in 1:time)
    {
        for (j in 1:time)
        {
            R[i,j]<-rho[abs(i-j)+1]
        }
    }
    return(R)
}

roo<-function(rho,time,corstr)
{
    if (rho>0.99)
    {rho<-0.99}
    roo<-rep(0,(time-1))
    if (corstr=="ar1")
    {
        for (i in 1:(time-1))
        {roo[i]<-rho^i}
    }
    else if (corstr=="exchangeable")
    {
        roo<-rep(rho,(time-1))
    }
    return(roo)
}




#'@title Calculate conditional probabilities for observing records at each time point
#'@description A function calculates conditional probabilities for longitudinal missing data. The observing probability is at observation-level.
#'@usage cond.prob(x_mis,gamma,id,time)
#'@param x_mis A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.
#'@param gamma coefficients calculated from missing data model
#'@param id A vector indicating subject id.
#'@param time The number of observations in total for each subject


#'@return a vector containing conditional probabilities.
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeesimdata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,scale = NULL,
#'          mismodel =obs_ind~x_mis1)
#'beta<-as.vector(summary(fit)$beta)
#'rho<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeesimdata$x_mis
#'pi<-cond.prob(x_mis,gamma,id,time=3)
#'
#'@export
cond.prob<-function(x_mis,gamma,id,time)
{
    samplesize<-length(unique(id))
    pi<-exp(x_mis%*%gamma)/(1+exp(x_mis%*%gamma))
    pi[rep(seq_len(time),time=samplesize)==1]<-1
    pi
}




#'@title Calculate the inverse of marginal probability for observing records at each time point
#'@description A function calculates the inverse of joint probabilities for weight calculation involved in WGEE.
#'@usage marg.prob(pi)
#'@param pi A matrix containing covariates for the missing data model. The first column should be all ones corresponding to the intercept.

#'@return a vector containing the inverse of joint probabilities.
#'@examples
#'## tests
#'# load data
#'data(wgeesimdata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeesimdata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,scale = NULL,
#'          mismodel =obs_ind~x_mis1)
#'beta<-as.vector(summary(fit)$beta)
#'rho<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeesimdata$x_mis
#'pi<-cond.prob(x_mis,gamma,id,time=3)
#'joint_prob<-marg.prob(pi)
#'@export
#'
marg.prob<-function(pi) #missing prob
{ pi[which(pi<0.00001)]<-0.00001
pii<-1/pi
time<-length(pi)
for (i in 2:time)
{
    pii[i]<-pii[i-1]*pii[i]
}
return(pii)
}

print.elcic.glm <- function(x, ...)
{
    cat ("ELCIC for glm\n")
    cat("****************\n")
    colnames(x) <- paste("model", 1:ncol(x))
    print.default(x[,], quote=FALSE, ...)
    cat("****************\n")
    for (i in 1:ncol(x))
        {cat("model ", i, ": ", deparse(attr(x, "formula")[[i]]), "\n", sep="")}
    cat("****************\n")
    cat("The model selected by ELCIC", ": ", deparse(attr(x, "formula")[[which.min(x[1,])]]), "\n", sep="")
    invisible()
}


print.elcic.gee <- function(x,y, ...)
{
    cat ("ELCIC for gee\n")
    cat("****************\n")
    colnames(x) <- paste("model", 1:ncol(x))
    print.default(x[,], quote=FALSE, ...)
    cat("****************\n")
    for (i in 1:ncol(x))
    {cat("model ", i, ": ", deparse(attr(x, "formula")[[i]]), "\n", sep="")}
    cat("****************\n")
    cat("The mean model selected by ELCIC", ": ",
        deparse(attr(x, "formula")[[which(x==min(x),arr.ind = TRUE)[2]]]), "\n", sep="")
    cat("The correlation structure selected by ELCIC",
        ": ", deparse(y[which(x==min(x),arr.ind = TRUE)[1]]), "\n", sep="")
    invisible()
}

print.elcic.wgee <- function(x,y, ...)
{
    cat ("ELCIC for wgee\n")
    cat("****************\n")
    colnames(x) <- paste("model", 1:ncol(x))
    print.default(x[,], quote=FALSE, ...)
    cat("****************\n")
    for (i in 1:ncol(x))
    {cat("model ", i, ": ", deparse(attr(x, "formula")[[i]]), "\n", sep="")}
    cat("****************\n")
    cat("The mean model selected by ELCIC", ": ",
        deparse(attr(x, "formula")[[which(x==min(x),arr.ind = TRUE)[2]]]), "\n", sep="")
    cat("The correlation structure selected by ELCIC",
        ": ", deparse(y[which(x==min(x),arr.ind = TRUE)[1]]), "\n", sep="")
    invisible()
}

print.mlic.wgee <- function(x,y, ...)
{
    cat ("MLIC for wgee\n")
    cat("****************\n")
    colnames(x) <- paste("model", 1:ncol(x))
    zz<-as.matrix(x[1,])
    colnames(zz)<-rownames(x)
    print.default(zz, quote=FALSE)
    cat("****************\n")
    for (i in 1:ncol(x))
    {cat("model ", i, ": ", deparse(attr(x, "formula")[[i]]), "\n", sep="")}
    cat("****************\n")
    cat("The mean model selected by MLIC", ": ",
        deparse(attr(x, "formula")[[which(x==min(x),arr.ind = TRUE)[2]]]), "\n", sep="")
    cat("The correlation structure selected by MLIC",
        ": ", deparse(rownames(x)), "\n", sep="")
    invisible()
}



print.qicw.wgee <- function(x,y, ...)
{
    cat ("QICW for wgee\n")
    cat("****************\n")
    colnames(x) <- paste("model", 1:ncol(x))
    zz<-as.matrix(x[1,])
    colnames(zz)<-rownames(x)
    print.default(zz, quote=FALSE)
    cat("****************\n")
    for (i in 1:ncol(x))
    {cat("model ", i, ": ", deparse(attr(x, "formula")[[i]]), "\n", sep="")}
    cat("****************\n")
    cat("The mean model selected by QICW", ": ",
        deparse(attr(x, "formula")[[which(x==min(x),arr.ind = TRUE)[2]]]), "\n", sep="")
    cat("The correlation structure selected by QICW",
        ": ", deparse(rownames(x)), "\n", sep="")
    invisible()
}


