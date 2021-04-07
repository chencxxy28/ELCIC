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

gee.generator<-function(beta,samplesize,time,num.time.dep,num.time.indep,rho=0.5,x.rho=0.2,dist,cor.str="exchangeable",x.cor.str="exchangeable") {
    if(cor.str=="exchangeable")
    {
        sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
        {
            for (i in 1:(time))
                sigma.dep[i,] <- rho
        }
    }else if(cor.str=="ar1"){
        sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
        {
            for (i in 1:(time))
                sigma.dep[i,] <- rho^(abs((1:(time))-i))
        }
    }else{
        sigma.dep <- matrix(0,(time),(time))
    }
    diag(sigma.dep)<-1

    if(x.cor.str=="exchangeable")
    {
        x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
        {
            for (i in 1:(time))
                x.sigma.dep[i,] <- x.rho
        }
    }else if(x.cor.str=="ar1"){
        x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
        {
            for (i in 1:(time))
                x.sigma.dep[i,] <- x.rho^(abs((1:(time))-i))
        }
    }else{
        x.sigma.dep <- matrix(0,(time),(time)) #covariance for the covariate matrix and residuals
    }
    diag(x.sigma.dep)<-1
    xf<-rep()
    y<-rep()

    #full X
    for (jj in 1:samplesize)
    {
    x.dep<-t(rmvnorm(num.time.dep,rep(0,time=time),x.sigma.dep))
    x.indep<-rep(rnorm(num.time.indep,0,1),time=time)
    xf.i<-cbind(1,x.dep,x.indep)
    #true X
    if(dist=="gaussian")
    {
        mu<-xf.i%*%beta
        y.i<-rmvnorm(1,mu,sigma.dep)
    }else if(dist=="poisson"){
        mu=exp(xf.i%*%beta)
        y.i<-genPoisNor(n=1,no.pois=length(mu),lamvec=mu,cmat.star=sigma.dep,no.norm=0,mean.vec=NULL,sd.vec=NULL)
    }else if(dist=="binomial"){
        mu=1/(1+exp(-xf.i%*%beta))
        y.i<-rmvbin(n=1,margprob=mu,bincorr=sigma.dep)
    }
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


R<-function(ro,id) #working correlation matrix
{
    time<-length(id)/length(unique(id))
    R<-diag(1,time,time)
    ro<-c(1,ro)
    for (i in 1:time)
    {
        for (j in 1:time)
        {
            R[i,j]<-ro[abs(i-j)+1]
        }
    }
    return(R)
}

roo<-function(ro,time,corstr)
{
    if (ro>0.99)
    {ro<-0.99}
    roo<-rep(0,(time-1))
    if (corstr=="ar1")
    {
        for (i in 1:(time-1))
        {roo[i]<-ro^i}
    }
    else if (corstr=="exchangeable")
    {
        roo<-rep(ro,(time-1))
    }
    return(roo)
}


