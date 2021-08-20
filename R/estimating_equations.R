

#'@title Estimating equation for ELCIC under GLM
#'@description A specified estimating equation for ELCIC under GLM. This estimating equation is used for marginal mean selection.
#'@usage ee.glm(x, y, betahat, dist)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param betahat A plug-in estimator solved by an external estimating procedure.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'
#'@return A matrix containing values of calculated estimating equations.
#'@examples
#'## tests
#'# load data
#'data(glmtoydata)
#'x<-glmtoydata$x
#'y<-glmtoydata$y
#'# obtain the estimates
#'fit<-glm(y~x-1,family="poisson")
#'betahat<-fit$coefficients
#'ee.matrix<-ee.glm(x, y, betahat, dist="poisson")
#'dim(ee.matrix)
#'
#'@export

######for cross-sectional glm
ee.glm<-function (x,y,betahat,dist)
{
    full<-ncol(x)
    samplesize<-nrow(x)
    ee<-matrix(0,nrow=full,ncol=samplesize)

    switch(dist,
           "gaussian"={    for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%betahat))
           }
           },
           "binomial"={ for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(1+exp(-t(as.matrix(x[i,],ncol=1))%*%betahat))^(-1))
           }
           },
           "poisson"={for (i in 1:samplesize)
           {
               ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%betahat))
           }
           },
           stop("Invalid type of dist. It should be one of gaussian,binomial,poisson")
    )

    ee

    # if(dist=="gaussian")
    # {
    # for (i in 1:samplesize)
    # {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(t(as.matrix(x[i,],ncol=1))%*%betahat))
    # }
    # }else if (dist=="poisson")
    # {
    #     for (i in 1:samplesize)
    #     {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-exp(t(as.matrix(x[i,],ncol=1))%*%betahat))
    #     }
    # }else if (dist=="binomial")
    # {
    #     for (i in 1:samplesize)
    #     {
    #     ee[,i]<-as.matrix(x[i,],ncol=1)%*%(y[i]-(1+exp(-t(as.matrix(x[i,],ncol=1))%*%betahat))^(-1))
    #     }
    # }

}




#######################################################
###### for joint selection of marginal mean and working correlation

#'@title Estimating equation for GEE without missingness or missing completely at random
#'@description Calculate estimating equation from GEE in ELCIC without missingness or missing completely at random. This estimating equation is used for joint selection of marginal mean and working correlation structure.
#'@usage ee.gee(y,x,r,id,beta,ro,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as GEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x<-geetoydata$x
#'y<-geetoydata$y
#'id<-geetoydata$id
#'corstr<-"exchangeable"
#'dist<-"poisson"
#'# obtain the estimates
#'library(geepack)
#'fit<-geeglm(y~x-1,data=geetoydata,family =dist,id=id,corstr = "ar1")
#'beta<-fit$coefficients
#'ro<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r<-rep(1,nrow(x))
#'ee.matrix<-ee.gee(y,x,r,id,beta,ro,phi,dist,corstr)
#'apply(ee.matrix,1,mean)
#'
#'@export
#'
######for longitudinal data gee
ee.gee<-function(y,x,r,id,beta,ro,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    time<-length(id)/n
    A<-diag(1,time,time)
    ro<-roo(ro,time,corstr)
    R<-R(ro,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    #z.col<-ncol(z)


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        #piii<-pii(pi[((i-1)*time+1):(i*time)])  #used when we have ipw
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        #WW<-W*r[((i-1)*time+1):(i*time)]*piii  #used when we have ipw
        WW<-W*r[((i-1)*time+1):(i*time)]
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        for (m in 1:(time-1)) #formula for ro
        {error<-0
        for (j in 1:(time-m))
        {
            error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        }
        error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        wgeei<-rbind(wgeei,error)
        }

        wgeef<-cbind(wgeef,wgeei)
    }

    return(wgeef)
}




#'@title Estimating equation for weighted GEE (WGEE) with data missing at random
#'@description Calculate estimating equation from WGEE in ELCIC with data missing at random. This estimating equation is used for joint selection of marginal mean and working correlation structure.
#'@usage ee.wgee(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes. use NA to indicate missing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records.
#'@param pi A vector containing observing probabilities across all observations.
#'@param time The number of observations for each subject.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as WGEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as WGEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeetoydata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,
#'      scale = NULL,mismodel =obs_ind~x_mis1+x_mis2)
#'beta<-as.vector(summary(fit)$beta)
#'ro<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeetoydata$x_mis
#'pi<-prob.obs(x_mis,gamma)
#'wgee.matrix<-ee.wgee(y=wgeetoydata$y,x=wgeetoydata$x,r=wgeetoydata$obs_ind,
#'pi=pi,id=wgeetoydata$id,time=3,beta=beta,ro=ro,phi=phi,dist=dist,corstr=corstr)
#'apply(wgee.matrix,1,mean)
#'
#'@export

ee.wgee<-function(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    ro<-roo(ro,time,corstr)
    A<-diag(1,time,time)
    R<-R(ro,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    # z.col<-ncol(z)
    # rlag<-ylag(id,r,1)
    # rlag[is.na(rlag)]<-1
    # S<-matrix(rep(0),nrow = z.col, ncol = n )


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        piii<-pii(pi[((i-1)*time+1):(i*time)])
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        WW<-W*r[((i-1)*time+1):(i*time)]*piii
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        for (m in 1:(time-1)) #formula for ro
        {error<-0
        for (j in 1:(time-m))
        {
            error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        }
        error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        wgeei<-rbind(wgeei,error)
        }

        wgeef<-cbind(wgeef,wgeei)
        # ss<-t(rlag[id==i]*(r[id==i]-lambda[id==i])*z[id==i,])
        # S[,i]<-as.matrix(apply(ss,1,sum))
    }
    return(wgeef)
}




#######################################################
###### only for marginal mean selection

#'@title Estimating equation of marginal mean for GEE without missing
#'@description Calculate estimating equation from GEE in ELCIC. This estimating equation is used for marginal mean selection.
#'@usage ee.gee.onlymean(y,x,r,id,beta,ro,phi,dist,corstr)
#'@param x A matrix containing covariates. The first column should be all ones corresponding to the intercept.
#'@param y A vector containing outcomes.
#'@param r A vector indicating missingness: 1 for observed records, and 0 for unobserved records. The default setup is that all data are observed.
#'@param id A vector indicating subject id.
#'@param beta A plug-in estimator solved by an external estimation procedure, such as GEE.
#'@param ro A correlation coefficients obtained from an external estimation procedure, such as GEE.
#'@param phi An over-dispersion parameter obtained from an external estimation procedure, such as GEE.
#'@param dist A specified distribution. It can be "gaussian", "poisson",and "binomial".
#'@param corstr A condidate correlation structure. It can be "independence","exchangeable", and "ar1".
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@note corstr should be prespecified.
#'
#'@examples
#'## tests
#'# load data
#'data(geetoydata)
#'x<-geetoydata$x
#'y<-geetoydata$y
#'id<-geetoydata$id
#'corstr<-"exchangeable"
#'dist<-"poisson"
#'# obtain the estimates
#'library(geepack)
#'fit<-geeglm(y~x-1,data=geetoydata,family =dist,id=id,corstr = corstr)
#'beta<-fit$coefficients
#'ro<-unlist(summary(fit)$corr[1])
#'phi<-unlist(summary(fit)$dispersion[1])
#'r<-rep(1,nrow(x))
#'ee.matrix<-ee.gee.onlymean(y,x,r,id,beta,ro,phi,dist,corstr)
#'apply(ee.matrix,1,mean)
#'
#'@export
ee.gee.onlymean<-function(y,x,r,id,beta,ro,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    time<-length(id)/n
    A<-diag(1,time,time)
    ro<-roo(ro,time,corstr)
    R<-R(ro,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    #z.col<-ncol(z)


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        #piii<-pii(pi[((i-1)*time+1):(i*time)])  #used when we have ipw
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        #WW<-W*r[((i-1)*time+1):(i*time)]*piii  #used when we have ipw
        WW<-W*r[((i-1)*time+1):(i*time)]
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        # e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        # for (m in 1:(time-1)) #formula for ro
        # {error<-0
        # for (j in 1:(time-m))
        # {
        #     error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        # }
        # error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        # wgeei<-rbind(wgeei,error)
        # }
        #
         wgeef<-cbind(wgeef,wgeei)
    }

    return(wgeef)
}




#'@title Estimating equation for marginal mean under weighted GEE(WGEE) with data missing at random
#'@description Calculate estimating equation from WGEE in ELCIC. This estimating equation is used for marginal mean selection.
#'@usage ee.wgee.onlymean(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
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
#'@note corstr should be prespecified.
#'
#'@return A matrix containing values of calculated estimating equations.
#'
#'@examples
#'## tests
#'# load data
#'data(wgeetoydata)
#'library(wgeesel)
#'data_wgee<-data.frame(do.call(cbind,wgeetoydata))
#'corstr<-"exchangeable"
#'dist<-"binomial"
#'id<-data_wgee$id
#'# obtain the estimates
#'fit<-wgee(y~x1+x2+x3,data_wgee,id,family=dist,corstr =corstr,
#'      scale = NULL,mismodel =obs_ind~x_mis1+x_mis2)
#'beta<-as.vector(summary(fit)$beta)
#'ro<-summary(fit)$corr
#'phi<-summary(fit)$phi
#'#calculate observing probabilies for all observations
#'gamma<-as.vector(summary(fit$mis_fit)$coefficients[,1])
#'x_mis<-wgeetoydata$x_mis
#'pi<-prob.obs(x_mis,gamma)
#'wgee.matrix<-ee.wgee.onlymean(y=wgeetoydata$y,x=wgeetoydata$x,r=wgeetoydata$obs_ind,
#'pi=pi,id=wgeetoydata$id,time=3,beta=beta,ro=ro,phi=phi,dist=dist,corstr=corstr)
#'apply(wgee.matrix,1,mean)
#'@export

ee.wgee.onlymean<-function(y,x,r,pi,id,time,beta,ro,phi,dist,corstr)
{ #full wgee
    n<-length(unique(id))
    p<-length(beta)
    ro<-roo(ro,time,corstr)
    A<-diag(1,time,time)
    R<-R(ro,id)
    W<-diag(1,time,time)
    V<-v(x,beta,dist)
    wgeef<-rep()

    # z.col<-ncol(z)
    # rlag<-ylag(id,r,1)
    # rlag[is.na(rlag)]<-1
    # S<-matrix(rep(0),nrow = z.col, ncol = n )


    y[which(is.na(y))]<-0
    for (i in 1:n)
    {
        piii<-pii(pi[((i-1)*time+1):(i*time)])
        AA<-A*V$v[((i-1)*time+1):(i*time)]^(-0.5)
        WW<-W*r[((i-1)*time+1):(i*time)]*piii
        e<-y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)]
        wgeei<-1/phi*t(V$der[((i-1)*time+1):(i*time),])%*%AA%*%solve(R)%*%AA%*%WW%*%e
        # e<-e*(V$v[((i-1)*time+1):(i*time)])^(-0.5)
        # for (m in 1:(time-1)) #formula for ro
        # {error<-0
        # for (j in 1:(time-m))
        # {
        #     error<-error+e[j]*e[j+m]*WW[j+m,j+m]
        # }
        # error<-error-ro[m]*phi*(n*(time-m)-p)/n #directly use phi might be more efficient
        # wgeei<-rbind(wgeei,error)
        # }

        wgeef<-cbind(wgeef,wgeei)
        # ss<-t(rlag[id==i]*(r[id==i]-lambda[id==i])*z[id==i,])
        # S[,i]<-as.matrix(apply(ss,1,sum))
    }
    return(wgeef)
}

