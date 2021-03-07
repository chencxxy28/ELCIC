##calculate I  #Now works for gaussian, poisson, and binomial distribution
II<-function (x,y,betahat,size,samplesize,dist)
{
    I_G<-matrix(0,nrow=size,ncol=size)
    if(dist=="gaussian")
    {
    for (i in 1:samplesize)
    {
        fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%betahat))
        I_G<-I_G+fitted%*%t(fitted)
    }
    }else if(dist=="poisson")
    {
        for (i in 1:samplesize)
        {
            fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%betahat))
            I_G<-I_G+fitted%*%t(fitted)
        }
    }else if(dist=="binomial")
    {
        for (i in 1:samplesize)
        {
            fitted<-as.matrix(x[i,],ncol=1)%*%(y[i]-1/(exp(-t(as.matrix(x[i,],ncol=1))%*%betahat)+1))
            I_G<-I_G+fitted%*%t(fitted)
        }
    }
    I_G
}

#calcuate J #Now works for gaussian, poisson, and binomial distribution
JJ<-function (x,y,betahat,size,samplesize,dist)
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
            fitted<-as.vector(exp(t(as.matrix(x[i,],ncol=1))%*%betahat))*as.matrix(x[i,],ncol=1)%*%t(as.matrix(x[i,],ncol=1))
            J_G<-J_G+fitted
        }
    }else if (dist=="binomial")
    {
        for (i in 1:samplesize)
        {
            pi_i<-1/(exp(-t(as.matrix(x[i,],ncol=1))%*%betahat)+1)
            fitted<-as.vector(pi*(1-pi))*as.matrix(x[i,],ncol=1)%*%t(as.matrix(x[i,],ncol=1))
            J_G<-J_G+fitted
        }
    }
    J_G
}

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
    if(dist=="gaussian")
    {
    y<-rnorm(samplesize,(x%*%beta),sd.gaussian)
    }else if(dist=="poisson")
    {
        y<-rpois(samplesize,exp(x%*%beta))
    }else if(dist=="binomial")
    {
        y<-rbinom(n=samplesize,size=1,prob=1/(1+exp(-x%*%beta)))
    }else if(dist=="NB")
    {
        y<-rnbinom(n=samplesize,size=ov,mu=exp(x%*%beta))
    }
    return(list(y=y,x=x))
}


