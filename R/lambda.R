##first derivative of -log EL
R1der<-function(lambda,ZZ)
{
    apply(ZZ,2,function(xx)
    {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

##second derivative of -log EL
R2der<-function(lambda,ZZ)
{
    r2der<-0
    for(i in 1:ncol(ZZ))
    {
        r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
        r2der<-r2der+r2der_i
    }
    r2der
}

##-log EL
R0der<-function(lambda,ZZ)
{
    apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

#function to find lambda, given beta and an
#'@title To calculate tuning parameter involved in ELCIC under GLM
#'@description This function aims to efficiently calculate the tuning parameter lambda in ELCIC.
#'@usage lambda.find.glm(x, y, betahat, dist)
#'@param x A matrix containing covariates. The first column should contain all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param betahat A plug-in estimator solved by an external estimating procedure.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'
#'@return A value of lambda (tuning parameter) vector involved in the empirical likelihood.
#'
#'@examples
#'## tests
#'# load data
#'data(glmtoydata)
#'x=glmtoydata$x
#'y=glmtoydata$y
#'# obtain the estimates
#'fit<-glm(y~x-1,family="poisson")
#'betahat<-fit$coefficients
#'lambda<-lambda.find.glm(x, y, betahat, dist="poisson")
#'lambda
#'
#'@export

lambda.find.glm<-function(x,y,betahat,dist)
{
    samplesize<-nrow(x)
    ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}





#'@title Calculate the tuning parameters involved in ELCIC under GEE
#'@description This function provides an efficient algorithm to calculate the tuning parameters involved in ELCIC under GEE.
#'@usage lambda.find.gee(x, y, id, betahat, r, dist, ro, phi, corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param betahat A plug-in estimator solved by an external estimation procedure, such as GEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return Tuning parameter values.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x=geetoydata$x
#'y=geetoydata$y
#'id=geetoydata$id
#'corstr="exchangeable"
#'dist="poisson"
#'# obtain the estimates
#'library(geepack)
#'fit<-geeglm(y~x-1,data=geetoydata,family =dist,id=id,corstr = corstr)
#'betahat<-fit$coefficients
#'ro<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r=rep(1,nrow(x))
#'lambda<-lambda.find.gee(x,y,id,betahat,r,dist,ro,phi,corstr)
#'lambda
#'
#'@export

#function to find lambda, given beta and an
lambda.find.gee<-function(x,y,id,betahat,r,dist,ro,phi,corstr)
{
    samplesize<-length(unique(id))

    ZZ<-ee.gee(y=y,x=x,r=r,id=id,beta=betahat,ro=ro,phi=phi,dist=dist,corstr=corstr)[,1:(samplesize)]
    #ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}




#'@title Calculate the tuning parameters involved in ELCIC under WGEE with data missing at random
#'@description This function provides an efficient algorithm to calculate the tuning parameters involved in ELCIC under WGEE with data missing at random.
#'@usage lambda.find.wgee(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param pi A vector containing observing probabilities across all observations.
#'@param time The number of observations for each subject
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as WGEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as WGEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return Tuning parameter values.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'library(wgeesel)
#'data_wgee=data.frame(do.call(cbind,wgeetoydata))
#'corstr="exchangeable"
#'dist="binomial"
#'id=data_wgee$id
#'# obtain the estimates
#'fit=wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,scale = NULL,
#'          mismodel =obs_ind~x_mis1+x_mis2)
#'beta=as.vector(summary(fit)$beta)
#'ro=summary(fit)$corr
#'phi=summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis=wgeetoydata$x_mis
#'pi=prob.obs(x_mis,gamma)
#'lambda=lambda.find.wgee(y=wgeetoydata$y,x=wgeetoydata$x,r=wgeetoydata$obs_ind,
#'pi=pi,id=wgeetoydata$id,time=3,beta=beta,ro=ro,phi=phi,dist=dist,corstr=corstr)
#'lambda
#'
#'@export

lambda.find.wgee<-function(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
{
    samplesize<-length(unique(id))

    ZZ<-ee.wgee(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)[,1:(samplesize)]
    #ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}





#function to find lambda under marginal mean selection in gee
#'@title Calculate the tuning parameters under marginal mean selection in gee
#'@description This function provides an efficient algorithm to calculate the tuning parameters involved in marginal mean selection in gee.
#'@usage lambda.find.gee.onlymean(x, y, id, betahat, r, dist, ro, phi, corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param betahat A plug-in estimator solved by an external estimation procedure, such as GEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return Tuning parameter values.
#'
#'@note corstr should be prespecified.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x=geetoydata$x
#'y=geetoydata$y
#'id=geetoydata$id
#'corstr="exchangeable"
#'dist="poisson"
#'# obtain the estimates
#'library(geepack)
#'fit<-geeglm(y~x-1,data=geetoydata,family =dist,id=id,corstr = corstr)
#'betahat<-fit$coefficients
#'ro<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r=rep(1,nrow(x))
#'lambda<-lambda.find.gee.onlymean(x,y,id,betahat,r,dist,ro,phi,corstr)
#'lambda
#'
#'@export
#'
lambda.find.gee.onlymean<-function(x,y,id,betahat,r,dist,ro,phi,corstr)
{
    samplesize<-length(unique(id))

    ZZ<-ee.gee.onlymean(y=y,x=x,r=r,id=id,beta=betahat,ro=ro,phi=phi,dist=dist,corstr=corstr)[,1:(samplesize)]
    #ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}






#'@title Calculate the tuning parameters involved in marginal mean selection under WGEE with data missing at random
#'@description This function provides an efficient algorithm to calculate the tuning parameters involved in marginal mean selection under WGEE with data missing at random.
#'@usage lambda.find.wgee.onlymean(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param pi A vector containing observing probabilities across all observations.
#'@param time The number of observations for each subject
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as WGEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as WGEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return Tuning parameter values.
#'
#'@note corstr should be prespecified.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'library(wgeesel)
#'data_wgee=data.frame(do.call(cbind,wgeetoydata))
#'corstr="exchangeable"
#'dist="binomial"
#'id=data_wgee$id
#'# obtain the estimates
#'fit=wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,scale = NULL,
#'          mismodel =obs_ind~x_mis1+x_mis2)
#'beta=as.vector(summary(fit)$beta)
#'ro=summary(fit)$corr
#'phi=summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis=wgeetoydata$x_mis
#'pi=prob.obs(x_mis,gamma)
#'lambda=lambda.find.wgee.onlymean(y=wgeetoydata$y,x=wgeetoydata$x,r=wgeetoydata$obs_ind,
#'pi=pi,id=wgeetoydata$id,time=3,beta=beta,ro=ro,phi=phi,dist=dist,corstr=corstr)
#'lambda
#'
#'@export

#function to find lambda under marginal mean selection in wgee
lambda.find.wgee.onlymean<-function(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
{
    samplesize<-length(unique(id))

    ZZ<-ee.wgee.onlymean(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)[,1:(samplesize)]
    #ZZ<-ee.glm(x=x,y=y,betahat=betahat,dist=dist)[,1:(samplesize)]
    #ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time)[,1:(n+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8

    repeat{
        rl<-R1der(lambda,ZZ)
        rll<-R2der(lambda,ZZ)
        Delta<--ginv(rll)%*%rl
        if(mean(abs(Delta))<tol)
        {break}else{
            repeat{
                mm<-0
                repeat{
                    delta<-gamma*Delta
                    index_1<-apply(ZZ,2,function (xx)
                    {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
                    )
                    if (sum(index_1)>0)
                    {gamma<-gamma/2
                    mm<-mm+1}else{break}}
                index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
                if (index_2==1)
                {gamma<-gamma/2}else{break}
            }
        }
        lambda<-lambda+delta
        c<-c+1
        gamma<-(c)^(-0.5)
    }
    lambda
}
